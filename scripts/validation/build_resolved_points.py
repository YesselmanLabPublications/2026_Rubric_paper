"""Per-molecule before/after values of the better-resolution cohort (Figure 4A-D).

One row per molecule of resolve_strict.csv, giving the value of every measure for the
worse-resolved (_w) and better-resolved (_b) structure. Suite, unassigned-suite and
pucker columns are taken from resolve_strict.csv. The pair-score columns are the mean
over the molecule's Review-tier pairs (score_w < 80) in resolve_pairs_matched.csv, and
n_review is their number (0 when the molecule has none).

Reads   data/reference/resolved_again/resolve_strict.csv        (build_resolved_cohort.py)
        data/reference/resolved_again/resolve_pairs_matched.csv (time_reversal_resolved.py)
Writes  data/reference/resolved_again/resolve_per_molecule_points.csv

Run after build_resolved_cohort.py and time_reversal_resolved.py.
"""
from __future__ import annotations

import argparse
import io
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _archive as A  # noqa: E402

RESOLVED = A.DATA / "reference" / "resolved_again"
COLUMNS = ["rep", "worst", "best", "gap", "method", "n_res", "n_review", "score_w", "score_b",
           "suit_w", "suit_b", "bang_w", "bang_b", "puck_n_w", "puck_n_b"]
KEY = ["rep", "worst", "best"]
REVIEW_CUTOFF = 80


def build(strict: Path, pairs: Path) -> pd.DataFrame:
    # Round, then pass through CSV text once, so the columns kept from resolve_strict.csv
    # carry the float representation of a written CSV.
    table = pd.read_csv(io.StringIO(pd.read_csv(strict)[COLUMNS].round(5).to_csv(index=False)))
    matched = pd.read_csv(pairs)
    review = (matched[matched.score_w < REVIEW_CUTOFF].groupby(KEY)
              .agg(score_w=("score_w", "mean"), score_b=("score_b", "mean"),
                   n_review=("score_w", "size")).reset_index())
    table = table.drop(columns=["score_w", "score_b", "n_review"]).merge(review, on=KEY, how="left")
    table["n_review"] = table.n_review.fillna(0).astype(int)
    return table[COLUMNS]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--strict", type=Path, default=RESOLVED / "resolve_strict.csv")
    ap.add_argument("--pairs", type=Path, default=RESOLVED / "resolve_pairs_matched.csv")
    ap.add_argument("--out", type=Path, default=RESOLVED / "resolve_per_molecule_points.csv")
    args = ap.parse_args()
    table = build(args.strict, args.pairs)
    table.to_csv(args.out, index=False)
    print(f"{len(table)} molecules, {(table.n_review > 0).sum()} with Review-tier pairs -> {args.out}")


if __name__ == "__main__":
    main()
