"""Unassignable-suite rate against resolution and wwPDB clashscore (Fig. S2B and C).

The structures are the entries of the curated and non-curated reference sets with a
resolution of 3.5 A or better and at least 20 suites (a rate from a handful of suites is
noise: one suite then moves a structure's rate by at most 5 points). Each structure's suites and
unassignable ('!!') suites are counted from the engine output in
data/reference/engine_archive/residues_orig.parquet. Within each resolution or
clashscore bin, structures are grouped into lineage clusters, so a molecule solved many
times counts once:

  rate         mean over clusters of each cluster's pooled rate (unassignable / suites, %)
  lo, hi       95% bootstrap interval of `rate` over clusters (4,000 resamples, seed 42)
  sd, sem      standard deviation over clusters of the mean per-structure rate, and sd/sqrt(n)
  pooled       unassignable / suites over all structures in the bin (weights by size;
               kept for comparison, not plotted)
  x            plotting position of the bin

The lineage clusters (data/reference/quality_axes/lineage_clusters.csv) group entries
whose base-pair fingerprints share a Jaccard index of 0.8 or more; an entry without a
cluster is its own cluster.

Reads   data/reference/engine_archive/residues_orig.parquet (build_engine_archive.py)
        data/reference/reference_sets/reference_set.json, reference_set_poor.json
        data/reference/quality_axes/lineage_clusters.csv
Writes  data/reference/quality_axes/bang_vs_quality.csv
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _archive as A  # noqa: E402

AXES = A.DATA / "reference" / "quality_axes"
SETS = A.DATA / "reference" / "reference_sets"
MAX_RESOLUTION = 3.5
MIN_SUITES = 20
RESAMPLES = 4000
SEED = 42

# axis: (column, [(lower bound exclusive, upper bound inclusive, label, x)])
BINS = {
    "resolution": ("resolution", [(-np.inf, 1.8, "<=1.8", 1.65), (1.8, 2.0, "1.8-2.0", 1.9),
                                  (2.0, 2.2, "2.0-2.2", 2.1), (2.2, 2.5, "2.2-2.5", 2.35),
                                  (2.5, 3.0, "2.5-3.0", 2.75), (3.0, 3.5, "3.0-3.5", 3.25)]),
    "clashscore": ("clashscore", [(-np.inf, 2, "<=2", 1.4), (2, 5, "2-5", 3.2), (5, 10, "5-10", 7.1),
                                  (10, 20, "10-20", 14.1), (20, np.inf, ">20", 27.0)]),
}


def structures(lineage: Path) -> pd.DataFrame:
    entries = {}
    for name in ("reference_set_poor.json", "reference_set.json"):
        entries.update(json.loads((SETS / name).read_text())["pdbs"])
    table = pd.DataFrame.from_dict(entries, orient="index")[["resolution", "clashscore"]]
    table = table.apply(pd.to_numeric, errors="coerce")
    residues = A.residues_orig()
    counts = residues.groupby("pdb_id").agg(n=("is_outlier", "size"), o=("is_outlier", "sum"))
    table = table.join(counts, how="inner")
    table = table[(table.resolution <= MAX_RESOLUTION) & (table.n >= MIN_SUITES)].copy()
    clusters = pd.read_csv(lineage).set_index("pdb_id").cluster
    ids = table.index.to_series()
    table["cluster"] = ids.map(clusters).fillna("u_" + ids)  # label order sets the bootstrap draws
    table["rate"] = 100 * table.o / table.n
    return table


def build(table: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)  # one generator, drawn bin after bin in table order
    rows = []
    for axis, (column, bins) in BINS.items():
        for lower, upper, label, x in bins:
            d = table[(table[column] > lower) & (table[column] <= upper)]
            by_cluster = d.groupby("cluster")
            pooled = (100 * by_cluster.o.sum() / by_cluster.n.sum()).to_numpy()
            boot = np.array([rng.choice(pooled, len(pooled)).mean() for _ in range(RESAMPLES)])
            lo, hi = np.percentile(boot, [2.5, 97.5])
            sd = by_cluster.rate.mean().std()
            rows.append({"axis": axis, "bin": label, "x": x, "rate": pooled.mean(), "lo": lo, "hi": hi,
                         "pooled": 100 * d.o.sum() / d.n.sum(), "n_molecules": len(pooled),
                         "n_struct": len(d), "n_suites": int(d.n.sum()), "sd": sd,
                         "sem": sd / np.sqrt(len(pooled))})
    return pd.DataFrame(rows).round(4)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lineage", type=Path, default=AXES / "lineage_clusters.csv")
    ap.add_argument("--out", type=Path, default=AXES / "bang_vs_quality.csv")
    args = ap.parse_args()
    table = structures(args.lineage)
    out = build(table)
    out.to_csv(args.out, index=False)
    print(f"{len(table):,} structures, {int(table.n.sum()):,} suites -> {args.out}")


if __name__ == "__main__":
    main()
