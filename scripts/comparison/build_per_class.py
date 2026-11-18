"""Per-class RUBRIC-DSSR comparison files, data/comparison/<class>.json (Fig. 2D/E, Fig. S3C/D).

Runs over every entry of the non-redundant set (uniqueRNAS.csv) that has DSSR output.
Inputs are the DSSR output in DSSR JSON layout (dssr_json.py), the deposited mmCIF files
and the RUBRIC engine: its normal output gives the selected pairs, `dump-classified` every
classified candidate and `dump-selection` the selection disposition of each candidate.
Afterwards run annotate_dssr_candidates.py (in this folder), which adds RUBRIC's class
for each DSSR pair RUBRIC did not report.
"""
import argparse
import csv
import json
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from functools import partial
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "data_gen"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _engine  # noqa: E402
import dssr_class_comparison as dcc  # noqa: E402

CLASSES = ("cWW", *dcc.NONCWW_CLASSES)


def strip_ambiguous(lw):
    if lw and lw.startswith("ambiguous (") and lw.endswith(")"):
        return lw[len("ambiguous ("):-1]
    return lw


def run_engine(binary, cif):
    """RUBRIC's selected pairs and classified candidates for one structure, or None."""
    try:
        out = subprocess.run([str(binary), str(cif), "--no-download"],
                             capture_output=True, text=True, timeout=600)
        data = json.loads(out.stdout)
        classified = subprocess.run([str(binary), "dump-classified", str(cif)],
                                    capture_output=True, text=True, timeout=600)
    except (subprocess.SubprocessError, ValueError):
        return None
    disposition = {}
    try:
        dump = subprocess.run([str(binary), "dump-selection", str(cif)],
                              capture_output=True, text=True, timeout=600)
        for line in dump.stdout.splitlines():
            t = line.split("\t")
            if len(t) >= 4:
                disposition[frozenset((t[0], t[1]))] = t[3]
    except subprocess.SubprocessError:
        pass
    selected = [{"r1": p["res_id1"], "r2": p["res_id2"],
                 "lw": strip_ambiguous(p["classification"]["lw_class"]),
                 "seq": p.get("sequence", ""), "cat": p.get("pair_category", "")}
                for p in data.get("pairs", [])]
    candidates = []
    for line in classified.stdout.splitlines():
        t = line.split("\t")
        if len(t) < 9:
            continue
        try:
            quality = float(t[8])
        except ValueError:
            quality = 0.0
        candidates.append(SimpleNamespace(
            res_id1=t[0], res_id2=t[1], lw_class=strip_ambiguous(t[2]),
            quality_score=quality, disposition=disposition.get(frozenset((t[0], t[1])))))
    return selected, candidates


def compare_entry(pdb_id, binary, structures, dssr_dir):
    engine = run_engine(binary, structures / f"{pdb_id}.cif")
    if engine is None:
        return None
    selected, candidates = engine
    dssr_json = dssr_dir / f"{pdb_id}.json"
    try:
        all_dssr = dcc.load_dssr_all(dssr_json)
    except (OSError, ValueError):
        all_dssr = {}
    return {lw: dcc.compare_structure(pdb_id, lw, dssr_json, selected, candidates, all_dssr)
            for lw in CLASSES}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--binary", type=Path, default=Path(_engine.default_engine()))
    ap.add_argument("--structures", type=Path, default=ROOT / "data" / "structures")
    ap.add_argument("--dssr-json", type=Path,
                    default=ROOT / "data" / "comparison" / "program_runs" / "dssr_json")
    ap.add_argument("--out", type=Path, default=ROOT / "data" / "comparison")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    pdbs = []
    with open(ROOT / "data" / "reference" / "reference_sets" / "uniqueRNAS.csv") as fh:
        for row in csv.DictReader(fh):
            p = (row.get("unique_pdb_id") or "").strip().upper()
            if p and p not in pdbs and (args.dssr_json / f"{p}.json").exists():
                pdbs.append(p)
    worker = partial(compare_entry, binary=args.binary.resolve(),
                     structures=args.structures.resolve(), dssr_dir=args.dssr_json.resolve())
    with ProcessPoolExecutor(args.workers) as ex:
        results = dict(zip(pdbs, ex.map(worker, pdbs, chunksize=4)))
    failed = [p for p, r in results.items() if r is None]

    args.out.mkdir(parents=True, exist_ok=True)
    for lw in CLASSES:
        per_entry = [results[p][lw] for p in sorted(results) if results[p] is not None]
        payload = {"summary": dcc.summarize(per_entry, lw),
                   "differences": [asdict(d) for r in per_entry for d in r.diffs]}
        (args.out / f"{lw}.json").write_text(json.dumps(payload, indent=2))
        s = payload["summary"]
        print(f"{lw}: {s['total_pdbs']} entries, DSSR pairs {s['total_dssr_pairs']}, "
              f"matched {s['matched']}")
    print(f"{len(pdbs)} entries, engine failures {len(failed)} {failed[:5]}")


if __name__ == "__main__":
    main()
