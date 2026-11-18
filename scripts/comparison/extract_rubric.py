"""RUBRIC's selected pairs, classified candidates and residue identities for every entry.

Uses the functions of extract_rubric_pairs.py, extract_rubric_candidates.py and
extract_residue_names.py (in this folder) on every entry of the
non-redundant set, and writes to data/comparison/program_runs/rubric/.
"""
import argparse
import csv
import hashlib
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import extract_residue_names as ern  # noqa: E402
import extract_rubric_candidates as erc  # noqa: E402
import extract_rubric_pairs as erp  # noqa: E402


def keyed(row):
    c1, n1, c2, n2, lw = row[:5]
    a, b = (c1, str(n1)), (c2, str(n2))
    labels = (lw or "?").split("|")
    if b < a:
        a, b = b, a
        labels = [x[0] + x[2] + x[1] if len(x) == 3 else x for x in labels]
    return (a, b), frozenset(labels)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=ROOT / "data" / "comparison" / "program_runs" / "rubric")
    ap.add_argument("--structures", type=Path, default=ROOT / "data" / "structures")
    ap.add_argument("--workers", type=int, default=6)
    args = ap.parse_args()
    binary = erp.DEFAULT_BINARY.resolve()
    with open(ROOT / "data" / "reference" / "reference_sets" / "uniqueRNAS.csv") as fh:
        pdbs = sorted({(r.get("unique_pdb_id") or "").strip().upper() for r in csv.DictReader(fh)} - {""})
    cif = lambda p: args.structures / f"{p}.cif"
    args.out.mkdir(parents=True, exist_ok=True)

    with ThreadPoolExecutor(args.workers) as ex:
        pairs = dict(zip(pdbs, ex.map(lambda p: erp.selected(binary, cif(p)), pdbs)))
    (args.out / "rubric_pairs.json").write_text(json.dumps(pairs))
    ok = sorted(p for p, r in pairs.items() if r["status"] == "ok")

    with ThreadPoolExecutor(args.workers) as ex:
        cands = dict(zip(ok, ex.map(lambda p: erc.candidates(binary, cif(p)), ok)))
    mismatched = [p for p in ok
                  if dict(keyed(r) for r in cands[p] if r[5] == "selected")
                  != dict(keyed(r) for r in pairs[p]["pairs"])]
    for p in mismatched:
        cands[p] = [[*r[:5], "selected"] for r in pairs[p]["pairs"]]
    data = {p: {"status": "ok", "pairs": cands[p]} for p in ok}
    data["_meta"] = {"fallback_to_output": mismatched,
                     "identified_dispositions": sorted(erc.IDENTIFIED),
                     "label_order": "input_residue_order",
                     "engine_sha256": hashlib.sha256(binary.read_bytes()).hexdigest()}
    (args.out / "rubric_candidates.json").write_text(json.dumps(data, separators=(",", ":")))

    with ThreadPoolExecutor(args.workers) as ex:
        names = dict(zip(pdbs, ex.map(lambda p: ern.residue_names(cif(p)), pdbs)))
    (args.out / "residue_names.json").write_text(json.dumps(names, separators=(",", ":")))
    print(f"{len(ok)} of {len(pdbs)} entries; selected pairs differing from dump-selection: {len(mismatched)}")


if __name__ == "__main__":
    main()
