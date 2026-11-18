"""Structure-level RUBRIC scores for every entry of the non-redundant set.

Runs the engine (`rubric --details`) on each entry of uniqueRNAS.csv and writes
data/reference/scores/scores_all_uniques.json, with one entry per structure:

    pairs_score          mean score of the scored pairs (0-100); null when the
                         structure has no scorable pair
    backbone_suiteness   mean suiteness of the scored suites (0-1)
    residues_score       backbone_suiteness x 100
    n_pairs, n_residues_scored, skipped_pairs, tier_summary, issue_summary, n_clashes

    python scripts/reference/score_structures_uniques.py --workers 6
"""

import argparse
import csv
import json
import subprocess
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts" / "data_gen"))
import _engine  # noqa: E402

DEFAULT_UNIQUES = REPO / "data" / "reference" / "reference_sets" / "uniqueRNAS.csv"
DEFAULT_OUT = REPO / "data" / "reference" / "scores" / "scores_all_uniques.json"


def read_ids(path):
    with open(path) as f:
        return [r["unique_pdb_id"].strip().upper() for r in csv.DictReader(f)
                if (r.get("unique_pdb_id") or "").strip()]


def score_one(pdb, engine):
    """(pdb, status, entry or error) for one structure."""
    args = [engine, *_engine.structure_args(pdb), "--details"]
    try:
        p = subprocess.run(args, capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired:
        return pdb, "timeout", None
    if p.returncode != 0 or not p.stdout.strip():
        return pdb, "failed", (p.stderr or "")[:200]
    try:
        d = json.loads(p.stdout)
    except json.JSONDecodeError:
        return pdb, "badjson", None

    pairs = d.get("pairs", [])
    scored = [p for p in pairs if isinstance(p.get("score"), (int, float))]
    issues = Counter(it["issue"] for pr in scored for it in (pr.get("issue_details") or [])
                     if it.get("severity", 0.0) > 0 and it.get("issue"))
    bb = d.get("backbone_suiteness")
    n_res = (d.get("score_details") or {}).get("n_residues_scored")
    if n_res is None:
        n_res = len(d.get("backbone_residues") or [])
    n_pairs = d.get("n_pairs") or 0
    entry = {
        "pairs_score": d.get("pairs_score") if n_pairs else None,
        "backbone_suiteness": bb,
        "residues_score": round(bb * 100.0, 2) if isinstance(bb, (int, float)) else None,
        "n_pairs": d.get("n_pairs"),
        "n_residues_scored": n_res,
        "skipped_pairs": len(pairs) - len(scored),
        "tier_summary": d.get("tier_summary") or {},
        "issue_summary": dict(issues),
        "n_clashes": d.get("n_clashes"),
    }
    return pdb, "ok", entry


def _stat(xs):
    s = sorted(x for x in xs if isinstance(x, (int, float)))
    if not s:
        return {}
    return {"n": len(s), "mean": round(sum(s) / len(s), 3),
            "median": round(s[len(s) // 2], 3), "min": round(s[0], 3), "max": round(s[-1], 3)}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--uniques", type=Path, default=DEFAULT_UNIQUES)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--engine", default=None, help="path to `rubric` (else RUBRIC_ENGINE / PATH)")
    ap.add_argument("--workers", type=int, default=6)
    args = ap.parse_args()

    engine = _engine.find_engine(args.engine)
    ids = read_ids(args.uniques)
    print(f"scoring {len(ids)} structures with {engine} ({args.workers} workers)", flush=True)

    scores, no_pairs, failed = {}, [], []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for pdb, status, payload in ex.map(partial(score_one, engine=engine), ids):
            if status == "ok":
                scores[pdb] = payload
                if not payload.get("n_pairs"):
                    no_pairs.append(pdb)
            else:
                failed.append({"pdb": pdb, "status": status, "error": payload})

    payload = {
        "scorer": "cpp:rubric",
        "rubric_bin": Path(engine).name,
        "summary": {
            "n_attempted": len(ids), "n_scored": len(scores),
            "n_with_pairs": sum(1 for e in scores.values()
                                if isinstance(e.get("pairs_score"), (int, float))),
            "n_no_pairs": len(no_pairs), "n_failed": len(failed),
            "pairs_score": _stat([e["pairs_score"] for e in scores.values()]),
            "backbone_suiteness": _stat([e["backbone_suiteness"] for e in scores.values()]),
        },
        "scores": scores,
        "no_pairs": sorted(set(no_pairs)),
        "failed": failed,
    }
    args.out.write_text(json.dumps(payload, indent=1))
    print(f"{len(scores)} scored, {len(no_pairs)} without pairs, {len(failed)} failed "
          f"({time.time() - t0:.0f} s) -> {args.out}")


if __name__ == "__main__":
    main()
