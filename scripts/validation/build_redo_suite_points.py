"""Suite and pucker columns of the PDB-REDO points table (Figure 4F-H).

For every structure in redo_per_structure_points.csv, the engine is run (`rubric
--details`) on the original coordinates and on the PDB-REDO re-refined coordinates.
Suites are matched between the two models by residue (chain, name and number), and
for the matched suites the script records

    n_res              number of matched suites
    suit_o, suit_r     mean suiteness, original and PDB-REDO
    bang_o, bang_r     unassigned suites (%), original and PDB-REDO
    puck_n_o, puck_n_r pucker outliers among the matched residues

the same per-molecule measures as the better-resolution cohort (build_resolved_cohort.py).
The pair-score columns are left as they are; build_review_points.py fills them.

Reads   data/reference/resolved_again/redo_per_structure_points.csv (structure list, clusters)
        data/validation/redo_files/<ID>/<ID>.pdb                   (PDB-REDO models)
        data/structures/<ID>.cif                                   (original models)
Writes  data/reference/resolved_again/redo_per_structure_points.csv

Run before build_review_points.py.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts" / "data_gen"))
import _engine  # noqa: E402

POINTS = REPO / "data" / "reference" / "resolved_again" / "redo_per_structure_points.csv"
REDO_FILES = REPO / "data" / "validation" / "redo_files"
SUITE_COLUMNS = ["n_res", "suit_o", "suit_r", "bang_o", "bang_r", "puck_n_o", "puck_n_r"]


def details(engine: str, args: list[str]) -> dict | None:
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "report.json"
        run = subprocess.run([engine, *args, "--details", "--out", str(out)],
                             capture_output=True, text=True)
        if run.returncode != 0 or not out.is_file():
            return None
        return json.loads(out.read_text())


def suites(report: dict) -> tuple[dict, set]:
    scored = {t["res_id"]: t for t in report.get("backbone_torsions") or []
              if t.get("scored") and t.get("suiteness") is not None}
    pucker = {p["res_id"] for p in report.get("pucker_outliers") or []}
    return scored, pucker


def structure_row(engine: str, pdb: str) -> dict:
    orig = details(engine, _engine.structure_args(pdb))
    redo = details(engine, [str(REDO_FILES / pdb / f"{pdb}.pdb")])
    if orig is None or redo is None:
        return {"pdb_id": pdb}
    (so, po), (sr, pr) = suites(orig), suites(redo)
    keys = sorted(set(so) & set(sr))
    n = len(keys)
    row = {"pdb_id": pdb, "n_res": n,
           "puck_n_o": sum(k in po for k in keys), "puck_n_r": sum(k in pr for k in keys)}
    if n:
        row.update(
            suit_o=sum(so[k]["suiteness"] for k in keys) / n,
            suit_r=sum(sr[k]["suiteness"] for k in keys) / n,
            bang_o=100 * sum(bool(so[k]["is_outlier"]) for k in keys) / n,
            bang_r=100 * sum(bool(sr[k]["is_outlier"]) for k in keys) / n,
        )
    return row


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--engine", default=None, help="path to `rubric` (else RUBRIC_ENGINE / PATH)")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--points", type=Path, default=POINTS)
    args = ap.parse_args()
    engine = _engine.find_engine(args.engine)

    points = pd.read_csv(args.points)
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        rows = list(ex.map(lambda p: structure_row(engine, p), points.pdb_id))
    new = pd.DataFrame(rows).set_index("pdb_id").reindex(columns=SUITE_COLUMNS)
    new[["suit_o", "suit_r", "bang_o", "bang_r"]] = new[["suit_o", "suit_r", "bang_o", "bang_r"]].round(5)
    columns = list(points.columns)
    points = points.drop(columns=SUITE_COLUMNS).join(new, on="pdb_id")
    for c in ("n_res", "puck_n_o", "puck_n_r"):
        points[c] = points[c].astype("Int64")
    points[columns].to_csv(args.points, index=False)
    failed = int(new.n_res.isna().sum())
    print(f"{len(points)} structures, {failed} without engine output -> {args.points}")


if __name__ == "__main__":
    main()
