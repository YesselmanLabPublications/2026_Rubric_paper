"""Engine output behind the better-resolution cohort and Fig. S2B/C (data/reference/engine_archive/).

Runs `rubric --details` once per structure and writes, for two sets of entries,

  pairs_<set>.parquet     one row per selected base pair: pdb_id, res_id1, res_id2,
                          lw_class, score, tier
  residues_<set>.parquet  one row per scored backbone suite: pdb_id, res_id, chain, num,
                          base_type, suiteness, is_outlier, pucker_outlier, pperp, delta

The residue rows are the suites that enter the engine's backbone suiteness (RNA
nucleotides whose C2' is bonded to an O2'; see the engine's has_rna_suite), so DNA and
residues modelled without their O2' are absent.

Sets
  orig    the non-redundant entries (unique_pdb_id of reference_sets/uniqueRNAS.csv)
  expand  further members of their BGSU equivalence classes, listed in
          engine_archive/expand_ids.txt (a fixed input)

Structures are read from data/structures/<ID>.cif; an entry that is not there is read
from, or downloaded into, data/structures/engine_cache/.
"""
import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "data_gen"))
import _engine  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
ARCHIVE = ROOT / "data" / "reference" / "engine_archive"
PAIR_COLUMNS = ["pdb_id", "res_id1", "res_id2", "lw_class", "score", "tier"]
RESIDUE_COLUMNS = ["pdb_id", "res_id", "chain", "num", "base_type", "suiteness",
                   "is_outlier", "pucker_outlier", "pperp", "delta"]
NUM = re.compile(r"-?\d+")


def run(engine, pdb):
    fd, out = tempfile.mkstemp(suffix=".json")
    os.close(fd)
    try:
        r = subprocess.run([engine, *_engine.structure_args(pdb), "--details", "--out", out],
                           capture_output=True, text=True)
        if r.returncode != 0 or os.path.getsize(out) == 0:
            return pdb, None
        return pdb, json.loads(Path(out).read_text())
    finally:
        Path(out).unlink(missing_ok=True)


def rows(pdb, d):
    pairs = [(pdb, p["res_id1"], p["res_id2"], (p.get("classification") or {}).get("lw_class"),
              p.get("score"), p.get("tier")) for p in d.get("pairs", [])]
    pucker = {o["res_id"]: o["pperp"] for o in d.get("pucker_outliers") or []}
    residues = []
    for t in d.get("backbone_torsions") or []:
        if not t.get("scored"):
            continue
        rid = t["res_id"]
        chain, num = rid.split("-")[0], NUM.match(rid.rsplit("-", 1)[-1])
        residues.append((pdb, rid, chain, int(num.group()) if num else None, t["base_type"],
                         t["suiteness"], t["is_outlier"], rid in pucker, pucker.get(rid),
                         t["delta"]))
    return pairs, residues


def build(engine, pdbs, workers):
    pairs, residues, failed = [], [], []
    with ThreadPoolExecutor(workers) as ex:
        for pdb, d in ex.map(lambda p: run(engine, p), pdbs):
            if d is None:
                failed.append(pdb)
                continue
            p, r = rows(pdb, d)
            pairs += p
            residues += r
    P = pd.DataFrame(pairs, columns=PAIR_COLUMNS)
    for c in ("pdb_id", "res_id1", "res_id2", "lw_class", "tier"):
        P[c] = P[c].astype("string")
    P["score"] = P["score"].astype("float64")
    R = pd.DataFrame(residues, columns=RESIDUE_COLUMNS)
    for c in ("pdb_id", "res_id", "chain", "base_type"):
        R[c] = R[c].astype("string")
    R["num"] = R["num"].astype("Int64")
    for c in ("suiteness", "pperp", "delta"):
        R[c] = R[c].astype("float64")
    for c in ("is_outlier", "pucker_outlier"):
        R[c] = R[c].astype(bool)
    return P, R, failed


def entries(part):
    if part == "orig":
        u = pd.read_csv(ROOT / "data/reference/reference_sets/uniqueRNAS.csv")
        return sorted(u.unique_pdb_id.str.upper().unique())
    return [s.strip() for s in (ARCHIVE / "expand_ids.txt").read_text().split() if s.strip()]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--engine", default=None)
    ap.add_argument("--out-dir", type=Path, default=ARCHIVE)
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()
    engine = _engine.find_engine(args.engine)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for part in ("orig", "expand"):
        pdbs = entries(part)
        if args.limit:
            pdbs = pdbs[: args.limit]
        P, R, failed = build(engine, pdbs, args.workers)
        P.to_parquet(args.out_dir / f"pairs_{part}.parquet", index=False)
        R.to_parquet(args.out_dir / f"residues_{part}.parquet", index=False)
        print(f"{part}: {len(P):,} pairs and {len(R):,} suites from {len(pdbs):,} entries; "
              f"{len(failed)} failed {failed[:10]}")


if __name__ == "__main__":
    main()
