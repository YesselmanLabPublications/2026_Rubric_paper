"""DSSR's base-pair tables in DSSR JSON layout, for build_per_class.py.

`x3dna-dssr -i=model.cif --pair-only --json` gives for each pair exactly the columns of the
text table written by the Table S4 command (index, nt1, nt2, bp, name, Saenger, LW, DSSR).
This writes the kept text output in that layout to
data/comparison/program_runs/dssr_json/<PDB>.json.
"""
import argparse
import gzip
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "data" / "comparison" / "program_runs"
FIELDS = ("index", "nt1", "nt2", "bp", "name", "Saenger", "LW", "DSSR")


def pair_table(text):
    rows, expected, inside = [], None, False
    for line in text.splitlines():
        m = re.match(r"^List of (\d+) base pairs?$", line.strip())
        if m:
            expected, inside = int(m.group(1)), True
            continue
        if not inside:
            continue
        tok = line.split()
        if tok and tok[0] == "nt1":
            continue
        if not tok or not tok[0].isdigit():
            if rows:
                break
            continue
        if len(tok) != len(FIELDS):
            raise ValueError(f"unexpected row {line!r}")
        row = dict(zip(FIELDS, tok))
        row["index"] = int(row["index"])
        rows.append(row)
    if expected is not None and len(rows) != expected:
        raise ValueError(f"{len(rows)} rows, header says {expected}")
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--runs", type=Path, default=RUNS)
    args = ap.parse_args()
    out = args.runs / "dssr_json"
    out.mkdir(exist_ok=True)
    files = sorted((args.runs / "raw" / "dssr").glob("*.dssr.txt.gz"))
    for f in files:
        with gzip.open(f, "rt", errors="replace") as fh:
            rows = pair_table(fh.read())
        (out / f"{f.name.split('.')[0]}.json").write_text(json.dumps({"num_pairs": len(rows), "pairs": rows}))
    print(f"wrote {len(files)} files to {out}")


if __name__ == "__main__":
    main()
