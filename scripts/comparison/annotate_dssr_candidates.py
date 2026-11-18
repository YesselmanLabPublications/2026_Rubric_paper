"""Add RUBRIC's own classification of each DSSR pair it did not report to the per-class
DSSR comparison files (data/comparison/cWW.json, tWW.json, ...).

Every dssr_only record gets `rubric_lw` (the LW class RUBRIC assigned the candidate,
None if it was never classified or never a candidate) and `rubric_disposition` (its
selection disposition). Figure 2D credits a DSSR pair whose class RUBRIC assigned,
whether or not selection kept it; a pair RUBRIC left unclassified is not credited.
"""
import argparse
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
COMPARISON = ROOT / "data" / "comparison"
CLASSES = ["cWW", "tWW", "cWH", "tWH", "cWS", "tWS", "cHH", "tHH", "cHS", "tHS", "cSS", "tSS"]
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "data_gen"))
import _engine  # noqa: E402

DEFAULT_BINARY = Path(_engine.default_engine())
DEFAULT_STRUCTURES = ROOT / "data" / "structures"


def dispositions(binary, path):
    run = subprocess.run([str(binary), "dump-selection", str(path)],
                         capture_output=True, text=True, check=True)
    out = {}
    for line in run.stdout.splitlines():
        a, b, lw, disposition = line.split("\t")
        if disposition != "invalid":
            out[frozenset((a, b))] = (None if lw == "None" else lw, disposition)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("binary", nargs="?", type=Path, default=DEFAULT_BINARY)
    ap.add_argument("structures", nargs="?", type=Path, default=DEFAULT_STRUCTURES)
    ap.add_argument("--comparison", type=Path, default=COMPARISON)
    args = ap.parse_args()
    binary, structures = args.binary.resolve(), args.structures.resolve()
    files = {lw: args.comparison / f"{lw}.json" for lw in CLASSES}
    data = {lw: json.loads(f.read_text()) for lw, f in files.items()}
    pdbs = sorted({d["pdb_id"] for payload in data.values()
                   for d in payload["differences"] if d["source"] == "dssr_only"})

    with ThreadPoolExecutor(max_workers=8) as pool:
        found = dict(zip(pdbs, pool.map(
            lambda p: dispositions(binary, structures / f"{p}.cif"), pdbs)))

    for lw, payload in data.items():
        for d in payload["differences"]:
            if d["source"] != "dssr_only":
                continue
            cls, disposition = found[d["pdb_id"]].get(
                frozenset((d["res_id1"], d["res_id2"])), (None, None))
            d["rubric_lw"], d["rubric_disposition"] = cls, disposition
        files[lw].write_text(json.dumps(payload, indent=2))
    print(f"annotated {len(files)} comparison files from {len(pdbs)} structures")


if __name__ == "__main__":
    main()
