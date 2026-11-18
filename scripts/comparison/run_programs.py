"""Run the six annotation programs on every entry of the non-redundant set and keep their output.

Each program runs with the version and command of Table S4 (scripts/benchmark/
benchmark_direct.py:command_for) on the deposited mmCIF file. RNAView reads only the first
character of a chain ID, so for entries with longer chain IDs it runs on a copy whose chains
are renamed to single characters; the names are restored when its output is converted.
Entries with more nucleic-acid chains than single characters are run in groups of chains
that together cover every pair of chains.

The native output is kept whenever a program writes it, also when the program then exits
with an error; convert_outputs.py decides which outputs are complete.

Output (resumable) in data/comparison/program_runs/:
  raw/<program>/<PDB>[.g<size>_<n>].<ext>.gz   native output
  status.jsonl                                one line per entry and program
  run_metadata.json                           versions and commands
"""
import argparse
import csv
import gzip
import json
import os
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "benchmark"))
import benchmark_all as common  # noqa: E402
import benchmark_direct as bd  # noqa: E402
from _annotators import (ALIAS_POOL, TooManyChains, _atom_site_scan, _build_alias,  # noqa: E402
                         _group_cif, alias_cif_chains, chain_groups)

PROGRAMS = ("dssr", "fr3d", "rnaview", "barnaba", "rnapolis", "bpnet")
# RNAView and BPNet fail on long file paths, so work directories are made in /tmp, as in
# the runtime benchmark.
WORK_DIR = "/tmp"
OUT = ROOT / "data" / "comparison" / "program_runs"


def native_output(program, cif, work, stdout):
    """(bytes, extension) of the program's base-pair output, or None if it wrote none."""
    if program == "dssr":
        return stdout.encode(), "dssr.txt"
    if program == "fr3d":
        for suffix in ("basepair", "basepairs"):
            p = work / f"{cif.stem}_{suffix}.txt"
            if p.exists():
                return p.read_bytes(), "fr3d_basepair.txt"
        return None
    path, ext = {"rnaview": (work / f"{cif.name}.out", "rnaview.out"),
                 "barnaba": (work / "annotation.ANNOTATE.pairing.out", "ANNOTATE.pairing.out"),
                 "rnapolis": (work / "annotation.json", "json"),
                 "bpnet": (work / f"{cif.stem}.out", "out")}[program]
    return (path.read_bytes(), ext) if path.exists() else None


def execute(program, cif, work, config, refdir, timeout):
    """Run one program on cif inside work; (return code or 'timeout', stdout, stderr, seconds)."""
    command = bd.command_for(program, cif, work, config)
    env = {**os.environ, **common.THREAD_ENV, "PYTHONHASHSEED": "0", "TMPDIR": str(work),
           "XDG_CACHE_HOME": str(work / "cache"), "MPLCONFIGDIR": str(work / "mpl"),
           "RNAVIEW": str(refdir / "rv"), "NUCLEIC_ACID_DIR": str(refdir / "bp") + "/"}
    t0 = time.perf_counter()
    proc = subprocess.Popen(command, cwd=work, env=env, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, start_new_session=True)
    try:
        out, err = proc.communicate(timeout=timeout)
        code = proc.returncode
    except subprocess.TimeoutExpired:
        common.terminate_tree(proc)
        out, err, code = b"", b"", "timeout"
    return code, out.decode(errors="replace"), err.decode(errors="replace"), time.perf_counter() - t0


def keep(out_dir, program, name, data, ext):
    dest = out_dir / "raw" / program / f"{name}.{ext}.gz"
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_suffix(".part")
    with gzip.open(part, "wb") as fh:
        fh.write(data)
    part.rename(dest)
    return str(dest.relative_to(out_dir))


def run_one(pdb, program, args, config, refdir):
    src = args.structures / f"{pdb}.cif"
    record = {"pdb": pdb, "program": program, "via": "direct", "returncodes": [], "outputs": []}
    t0 = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix="rubric-cmp-", dir=WORK_DIR) as temp:
        work = Path(temp)
        jobs = [(pdb, None)]
        if program == "rnaview":
            text = src.read_text(errors="replace")
            try:
                aliased = alias_cif_chains(text)
                if aliased is not None:
                    record["via"] = "renamed"
                    jobs = [(pdb, aliased[0])]
            except TooManyChains:
                return run_rnaview_groups(pdb, text, record, args, config, refdir, t0)
        for name, text in jobs:
            cif = work / f"{name}.cif"
            if text is None:
                cif.symlink_to(src)
            else:
                cif.write_text(text)
            code, stdout, stderr, _ = execute(program, cif, work, config, refdir, args.timeout)
            record["returncodes"].append(code)
            native = native_output(program, cif, work, stdout)
            if native is not None:
                record["outputs"].append(keep(args.out, program, name, *native))
            if code != 0:
                record["stderr_tail"] = stderr[-400:]
    record["seconds"] = round(time.perf_counter() - t0, 3)
    return record


def run_rnaview_groups(pdb, text, record, args, config, refdir, t0):
    """RNAView on groups of chains, halving the group size until every group completes."""
    lines = text.splitlines()
    auth_i, _atom_i, _chains, nucleic = _atom_site_scan(lines)
    size = len(ALIAS_POOL)
    while size >= 2:
        outputs, codes, complete = [], [], True
        with tempfile.TemporaryDirectory(prefix="rubric-cmp-", dir=WORK_DIR) as temp:
            work = Path(temp)
            for g, group in enumerate(chain_groups(sorted(nucleic), size)):
                cif = work / f"{pdb}_s{size}_g{g}.cif"
                cif.write_text(_group_cif(lines, auth_i, _build_alias(group, set(group))))
                code, stdout, _err, _ = execute("rnaview", cif, work, config, refdir, args.timeout)
                codes.append(code)
                native = native_output("rnaview", cif, work, stdout)
                if native is None or b"END_base-pair" not in native[0]:
                    complete = False
                    break
                outputs.append((f"{pdb}.g{size}_{g}", native))
            if complete:
                record.update(via=f"chain_groups_{size}", returncodes=codes,
                              outputs=[keep(args.out, "rnaview", n, *nat) for n, nat in outputs],
                              seconds=round(time.perf_counter() - t0, 3))
                return record
        size //= 2
    record.update(via="chain_groups", returncodes=codes, outputs=[],
                  seconds=round(time.perf_counter() - t0, 3))
    return record


def default_pdbs():
    with open(ROOT / "data" / "reference" / "reference_sets" / "uniqueRNAS.csv") as fh:
        ids = [(row.get("unique_pdb_id") or "").strip().upper() for row in csv.DictReader(fh)]
    return sorted({p for p in ids if p})


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pdbs", type=Path, help="file with one PDB ID per line (default: uniqueRNAS.csv)")
    ap.add_argument("--programs", default=",".join(PROGRAMS))
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--structures", type=Path, default=ROOT / "data" / "structures")
    ap.add_argument("--config", type=Path, default=ROOT / ".runtime-tools" / "config.json")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--timeout", type=float, default=1800)
    args = ap.parse_args()
    config = json.loads(args.config.read_text())
    pdbs = ([l.strip() for l in open(args.pdbs) if l.strip()] if args.pdbs else default_pdbs())
    programs = args.programs.split(",")
    args.out.mkdir(parents=True, exist_ok=True)
    status = args.out / "status.jsonl"
    done = set()
    if status.exists():
        done = {(r["pdb"], r["program"]) for r in map(json.loads, status.read_text().splitlines())}
    meta = {"programs": programs, "n_entries": len(pdbs),
            "commands_from": "scripts/benchmark/benchmark_direct.py:command_for",
            "packages": config.get("packages"), "git_commits": config.get("git_commits")}
    (args.out / "run_metadata.json").write_text(json.dumps(meta, indent=1))
    jobs = [(p, t) for p in pdbs for t in programs if (p, t) not in done]
    print(f"{len(jobs)} runs to do ({len(done)} already recorded)", flush=True)
    with tempfile.TemporaryDirectory(prefix="rubric-cmp-ref-", dir=WORK_DIR) as ref, \
            status.open("a") as log, ThreadPoolExecutor(args.workers) as ex:
        refdir = Path(ref)
        (refdir / "rv").symlink_to(config["rnaview_home"], target_is_directory=True)
        (refdir / "bp").symlink_to(config["bpnet_sysfiles"], target_is_directory=True)
        futures = [ex.submit(run_one, p, t, args, config, refdir) for p, t in jobs]
        for n, f in enumerate(as_completed(futures), 1):
            record = f.result()
            log.write(json.dumps(record) + "\n")
            log.flush()
            if n % 100 == 0:
                print(f"[{n}/{len(jobs)}]", flush=True)


if __name__ == "__main__":
    main()
