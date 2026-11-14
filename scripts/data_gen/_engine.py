"""Run the RUBRIC C++ engine (`rubric --details`) and return its JSON."""

import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
STRUCTURES = REPO / "data" / "structures"
ENGINE_CACHE = STRUCTURES / "engine_cache"


def structure_args(pdb_id: str) -> list[str]:
    local = STRUCTURES / f"{pdb_id.upper()}.cif"
    if local.exists():
        return [str(local), "--no-download"]
    return [pdb_id, "--cache-dir", str(ENGINE_CACHE)]


def _installed_binary() -> str | None:
    """The program bundled with the installed rubric-rna package, if any.

    Calling it directly skips the small Python launcher that `pip install` puts on
    PATH; the results are the same.
    """
    try:
        import rubric_rna
    except ImportError:
        return None
    binary = rubric_rna.binary_path()
    return str(binary) if binary.is_file() else None


def candidates(explicit: str | None = None) -> list[str | None]:
    """Where to look for the engine, in order: an explicit path, $RUBRIC_ENGINE,
    the installed rubric-rna package, then `rubric` on PATH."""
    return [explicit, os.environ.get("RUBRIC_ENGINE"), _installed_binary(), shutil.which("rubric")]


def find_engine(explicit: str | None = None) -> str:
    for cand in candidates(explicit):
        if cand and Path(cand).exists():
            return str(Path(cand).resolve())
    raise FileNotFoundError(
        "RUBRIC engine not found. Install it with `pip install -r requirements-engine.txt`, "
        "or set RUBRIC_ENGINE=/path/to/rubric."
    )


def default_engine() -> str:
    """find_engine() for use as a command-line default: never raises."""
    try:
        return find_engine()
    except FileNotFoundError:
        return "rubric"


def run_rubric(pdb_id: str, engine: str, cache_dir: str | None = None) -> dict:
    fd, out = tempfile.mkstemp(suffix=".json")
    os.close(fd)
    try:
        src = [pdb_id, "--cache-dir", str(cache_dir)] if cache_dir else structure_args(pdb_id)
        cmd = [engine, *src, "--details", "--out", out]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0 or os.path.getsize(out) == 0:
            raise RuntimeError(f"rubric failed for {pdb_id}: {r.stderr.strip()[:200]}")
        with open(out) as f:
            return json.load(f)
    finally:
        try:
            os.unlink(out)
        except OSError:
            pass


_DNA_LITERAL = frozenset({"DA", "DC", "DG", "DT", "DU"})


def is_dna_name(name: str) -> bool:
    n = (name or "").strip().upper()
    return n in _DNA_LITERAL


def is_dna_res_id(res_id: str) -> bool:
    return is_dna_name(res_name(res_id))


def is_scorable(pair: dict) -> bool:
    c = pair.get("classification") or {}
    lw = c.get("lw_class") or ""
    if pair.get("pair_category") != "base-base" or c.get("is_ambiguous"):
        return False
    if not lw or "|" in lw:
        return False
    if is_dna_res_id(pair.get("res_id1", "")) or is_dna_res_id(pair.get("res_id2", "")):
        return False
    if "score" in pair and pair["score"] is None:
        return False
    return True


def res_name(res_id: str) -> str:
    return res_id.split("-")[1]
