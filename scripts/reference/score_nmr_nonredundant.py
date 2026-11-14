"""RUBRIC scores for a non-redundant set of solution-NMR RNA structures (Figure S7).

NMR entries are not in uniqueRNAS.csv, so they get their own non-redundant set built
the same way: one entry per BGSU equivalence class. Candidates are solution-NMR PDB
entries with an RNA polymer and no protein or DNA polymer (RCSB search API); in each
BGSU class (all-resolution list) the first NMR candidate in BGSU's member order is
kept. The engine reads the first model, so each entry is scored on model 1.

The BGSU list and the chosen IDs are saved under data/reference/reference_sets/ and
reused on later runs; --refresh rebuilds them (RCSB and BGSU change over time).
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import tempfile
import urllib.request
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SETS = REPO / "data" / "reference" / "reference_sets"
STRUCTURES = REPO / "data" / "structures" / "nmr"
OUT = REPO / "data" / "reference" / "scores" / "scores_nmr_nonredundant.json"
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "data_gen"))
import _engine  # noqa: E402

ENGINE = _engine.default_engine()

BGSU_RELEASE = "4.62"
BGSU_URL = "https://rna.bgsu.edu/rna3dhub/nrlist/download/{release}/all/csv"
RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_FILE = "https://files.rcsb.org/download/{pdb}.cif"


def _terminal(attribute: str, operator: str, value) -> dict:
    return {"type": "terminal", "service": "text",
            "parameters": {"attribute": attribute, "operator": operator, "value": value}}


def nmr_rna_entries() -> set[str]:
    query = {"query": {"type": "group", "logical_operator": "and", "nodes": [
        _terminal("exptl.method", "exact_match", "SOLUTION NMR"),
        _terminal("rcsb_entry_info.polymer_entity_count_RNA", "greater_or_equal", 1),
        _terminal("rcsb_entry_info.polymer_entity_count_protein", "equals", 0),
        _terminal("rcsb_entry_info.polymer_entity_count_DNA", "equals", 0)]},
        "return_type": "entry", "request_options": {"return_all_hits": True}}
    req = urllib.request.Request(RCSB_SEARCH, data=json.dumps(query).encode(),
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req) as resp:
        return {hit["identifier"].upper() for hit in json.load(resp)["result_set"]}


def bgsu_list(release: str, refresh: bool) -> Path:
    path = SETS / f"bgsu_nrlist_all_{release}.csv"
    if refresh or not path.exists():
        urllib.request.urlretrieve(BGSU_URL.format(release=release), path)
    return path


def nonredundant(release: str, refresh: bool) -> list[tuple[str, str]]:
    """(pdb_id, bgsu_class) per class, saved to nmr_nonredundant.csv."""
    path = SETS / "nmr_nonredundant.csv"
    if path.exists() and not refresh:
        with open(path, newline="") as handle:
            return [(r["pdb_id"], r["bgsu_class"]) for r in csv.DictReader(handle)]
    nmr = nmr_rna_entries()
    chosen = {}
    with open(bgsu_list(release, refresh), newline="") as handle:
        for cls, _rep, members in csv.reader(handle):
            ids = [m.split("|")[0].upper() for m in members.split(",")]
            first = next((p for p in ids if p in nmr), None)
            if first and first not in chosen:
                chosen[first] = cls
    rows = sorted(chosen.items())
    with open(path, "w", newline="") as handle:
        w = csv.writer(handle)
        w.writerow(["pdb_id", "bgsu_class"])
        w.writerows(rows)
    print(f"{len(nmr)} NMR RNA entries -> {len(rows)} non-redundant (BGSU {release})")
    return rows


def fetch(pdb: str) -> None:
    cif = STRUCTURES / f"{pdb}.cif"
    if not cif.exists():
        urllib.request.urlretrieve(RCSB_FILE.format(pdb=pdb), cif)


def score_one(pdb: str) -> tuple[str, str, dict]:
    cif = STRUCTURES / f"{pdb}.cif"
    if not cif.exists():
        return pdb, "missing", {}
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as fh:
        out = Path(fh.name)
    try:
        run = subprocess.run([ENGINE, str(cif), "--no-download", "--out", str(out)],
                             capture_output=True, text=True, timeout=600)
        if run.returncode != 0:
            return pdb, "failed", {}
        d = json.loads(out.read_text())
    except (subprocess.TimeoutExpired, json.JSONDecodeError):
        return pdb, "failed", {}
    finally:
        out.unlink(missing_ok=True)
    return pdb, "ok", {
        "pairs_score": d.get("pairs_score"),
        "backbone_suiteness": d.get("backbone_suiteness"),
        "n_pairs": d.get("n_pairs", 0),
        "n_scored": sum(p.get("score") is not None for p in d.get("pairs", [])),
        "n_clashes": d.get("n_clashes", 0),
        "tier_summary": d.get("tier_summary"),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--release", default=BGSU_RELEASE, help="BGSU non-redundant list release")
    ap.add_argument("--refresh", action="store_true", help="rebuild the BGSU / RCSB lists")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    STRUCTURES.mkdir(parents=True, exist_ok=True)
    rows = nonredundant(args.release, args.refresh)
    classes = dict(rows)
    ids = [pdb for pdb, _ in rows]
    with ProcessPoolExecutor(args.workers) as ex:
        list(ex.map(fetch, ids))
        results = list(ex.map(score_one, ids, chunksize=4))

    status, scores = {}, {}
    for pdb, st, info in results:
        status[st] = status.get(st, 0) + 1
        if st == "ok":
            scores[pdb] = {**info, "bgsu_class": classes[pdb]}
    doc = {"scorer": "cpp:rubric", "bgsu_release": args.release,
           "selection": "solution NMR, RNA polymer, no protein/DNA; one entry per BGSU class; model 1",
           "summary": {"n_attempted": len(ids), **{f"n_{k}": v for k, v in status.items()}},
           "scores": scores}
    tmp = args.out.with_suffix(".json.part")
    tmp.write_text(json.dumps(doc, indent=1))
    tmp.replace(args.out)
    print(f"{len(ids):,} entries {status} -> {args.out.relative_to(REPO)}")


if __name__ == "__main__":
    main()
