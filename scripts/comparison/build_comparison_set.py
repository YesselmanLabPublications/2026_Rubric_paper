"""Comparison set: the entries that RUBRIC and all six programs processed.

Reads data/comparison/program_runs/pairs/ (convert_outputs.py) and
data/comparison/program_runs/rubric/ (extract_rubric.py) and writes, for the entries in which
every program has usable output, the files read by the figure notebooks and the validation
scripts:
  data/comparison/<program>_pairs.json, rubric_pairs.json, rubric_candidates.json,
  residue_names.json, comparison_entries.json and fr3d/<PDB>.csv (FR3D in the
  BGSU RNA 3D Hub CSV layout, "unit1","class","unit2").
"""
import argparse
import gzip
import json
import shutil
from pathlib import Path

from convert_outputs import fr3d_unit_in_model

ROOT = Path(__file__).resolve().parents[2]
COMPARISON = ROOT / "data" / "comparison"
PROGRAMS = ("dssr", "rnaview", "barnaba", "rnapolis", "bpnet", "fr3d_local")


def ok_entries(d):
    return {p for p, v in d.items() if v.get("status") == "ok"}


def fr3d_csv_rows(raw_gz):
    rows = []
    with gzip.open(raw_gz, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or not f[0]:
                continue
            if not (fr3d_unit_in_model(f[0]) and fr3d_unit_in_model(f[2])):
                continue
            residue = lambda u: (u[2], u[4], u[7] if len(u) >= 8 else "") if len(u) >= 5 else u
            if residue(f[0].split("|")) != residue(f[2].split("|")):
                rows.append(f'"{f[0]}","{f[1]}","{f[2]}"')
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--runs", type=Path, default=COMPARISON / "program_runs")
    ap.add_argument("--out", type=Path, default=COMPARISON)
    args = ap.parse_args()
    pairs = {p: json.loads((args.runs / "pairs" / f"{p}_pairs.json").read_text()) for p in PROGRAMS}
    rubric_dir = args.runs / "rubric"
    rubric = json.loads((rubric_dir / "rubric_pairs.json").read_text())
    # A two-class label means ambiguous; two DNA pairs in the RUBRIC output (7VBC, 8EHQ)
    # carry one without the flag.
    for v in rubric.values():
        for r in v["pairs"]:
            if "|" in (r[4] or "") and not r[5]:
                r[5] = True

    common = ok_entries(rubric)
    for d in pairs.values():
        common &= ok_entries(d)
    common = sorted(common)
    keep = lambda d: {p: d[p] for p in common}

    args.out.mkdir(parents=True, exist_ok=True)
    for program, d in pairs.items():
        (args.out / f"{program}_pairs.json").write_text(json.dumps(keep(d)))
    (args.out / "rubric_pairs.json").write_text(json.dumps(keep(rubric)))
    cands = json.loads((rubric_dir / "rubric_candidates.json").read_text())
    meta = cands.pop("_meta")
    meta["fallback_to_output"] = [p for p in meta["fallback_to_output"] if p in set(common)]
    (args.out / "rubric_candidates.json").write_text(
        json.dumps({**keep(cands), "_meta": meta}, separators=(",", ":")))
    names = json.loads((rubric_dir / "residue_names.json").read_text())
    (args.out / "residue_names.json").write_text(json.dumps(keep(names), separators=(",", ":")))
    (args.out / "comparison_entries.json").write_text(json.dumps(common))

    fr3d_dir = args.out / "fr3d"
    if fr3d_dir.exists():
        shutil.rmtree(fr3d_dir)
    fr3d_dir.mkdir()
    for p in common:
        rows = fr3d_csv_rows(args.runs / "raw" / "fr3d" / f"{p}.fr3d_basepair.txt.gz")
        (fr3d_dir / f"{p}.csv").write_text("\n".join(rows) + ("\n" if rows else ""))
    print(f"comparison set: {len(common)} entries")


if __name__ == "__main__":
    main()
