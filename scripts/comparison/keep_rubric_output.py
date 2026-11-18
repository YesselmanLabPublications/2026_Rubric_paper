"""RUBRIC's own output for every entry of the non-redundant set, kept beside the native output
of the six programs in data/comparison/program_runs/raw/rubric/<PDB>.rubric.json.gz.

Each file is the JSON written by `rubric <PDB>.cif --no-download`, the command whose output
extract_rubric.py reads. Entries already kept are skipped.
"""
import argparse
import csv
import gzip
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "data_gen"))
import _engine  # noqa: E402


def keep(engine, cif, dest):
    if dest.exists():
        return True
    run = subprocess.run([engine, str(cif), "--no-download"], capture_output=True)
    if run.returncode != 0:
        return False
    part = dest.with_suffix(".part")
    part.write_bytes(gzip.compress(run.stdout, mtime=0))
    part.rename(dest)
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--engine", default=None, help="path to `rubric` (else RUBRIC_ENGINE / PATH)")
    ap.add_argument("--out", type=Path, default=ROOT / "data" / "comparison" / "program_runs" / "raw" / "rubric")
    ap.add_argument("--structures", type=Path, default=ROOT / "data" / "structures")
    ap.add_argument("--workers", type=int, default=6)
    args = ap.parse_args()
    engine = _engine.find_engine(args.engine)
    with open(ROOT / "data" / "reference" / "reference_sets" / "uniqueRNAS.csv") as fh:
        pdbs = sorted({(r.get("unique_pdb_id") or "").strip().upper() for r in csv.DictReader(fh)} - {""})
    args.out.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(args.workers) as ex:
        ok = list(ex.map(lambda p: keep(engine, args.structures / f"{p}.cif",
                                        args.out / f"{p}.rubric.json.gz"), pdbs))
    failed = [p for p, k in zip(pdbs, ok) if not k]
    print(f"{len(pdbs) - len(failed)} of {len(pdbs)} entries kept in {args.out}"
          + (f"; failed: {' '.join(failed)}" if failed else ""))


if __name__ == "__main__":
    main()
