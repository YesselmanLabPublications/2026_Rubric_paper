"""Helpers shared by the runtime benchmark (benchmark_direct.py) and run_programs.py: structure
selection, trial records, process cleanup and the timing summaries."""
from collections import defaultdict
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import signal
import statistics
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ("rubric", "rubric_noscore", "dssr", "fr3d", "rnaview", "barnaba", "rnapolis", "bpnet")
THREAD_ENV = {key: "1" for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS")}


def sha(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def write_json(path, data):
    temp = path.with_suffix(path.suffix + ".part")
    temp.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    temp.replace(path)


def choose_structures(rows, limit, seed):
    if not limit or limit >= len(rows):
        return sorted(rows, key=lambda r: r["pdb"])
    if limit < 1:
        raise ValueError("limit must be nonnegative")
    buckets = [[] for _ in range(4)]
    for row in rows:
        size = row["nucleotides"]
        buckets[0 if size <= 50 else 1 if size <= 200 else 2 if size <= 1000 else 3].append(row)
    rng = random.Random(seed)
    for bucket in buckets:
        rng.shuffle(bucket)
    chosen = []
    while len(chosen) < limit:
        for bucket in buckets:
            if bucket and len(chosen) < limit:
                chosen.append(bucket.pop())
    return sorted(chosen, key=lambda r: r["pdb"])


def terminate_tree(proc):
    try:
        os.killpg(proc.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    proc.communicate()


def read_trials(path):
    if not path.exists():
        return []
    raw = path.read_bytes()
    if raw and not raw.endswith(b"\n"):
        end = raw.rfind(b"\n") + 1
        path.with_suffix(".interrupted-tail").write_bytes(raw[end:])
        path.write_bytes(raw[:end])
        raw = raw[:end]
    return [json.loads(line) for line in raw.splitlines() if line.strip()]


def hardware():
    info = {"cpu_count": os.cpu_count()}
    if sys.platform == "darwin":
        for name, key in (("cpu_model", "machdep.cpu.brand_string"), ("memory_bytes", "hw.memsize")):
            result = subprocess.run(["sysctl", "-n", key], capture_output=True, text=True)
            info[name] = result.stdout.strip() if result.returncode == 0 else "unavailable"
    elif Path("/proc/cpuinfo").exists():
        info["cpu_model"] = next((s.split(":", 1)[1].strip() for s in
            Path("/proc/cpuinfo").read_text().splitlines() if s.startswith("model name")), platform.machine())
        info["memory"] = Path("/proc/meminfo").read_text().splitlines()[0]
    return info


def summarize(trials, manifest, tools, repeats):
    grouped = defaultdict(list)
    for trial in trials:
        if not trial["warmup"]:
            grouped[(trial["pdb"], trial["tool"])].append(trial)
    summaries = []
    for row in manifest:
        for tool in tools:
            runs = grouped[(row["pdb"], tool)]
            good = [r for r in runs if r["status"] == "ok"]
            complete = len(runs) == repeats and len(good) == repeats
            signatures = {r.get("output_sha256") for r in good}
            stable = len(signatures) == 1
            eligible = complete and stable
            summaries.append({"pdb": row["pdb"], "nucleotides": row["nucleotides"], "tool": tool,
                              "attempted": len(runs), "successful": len(good),
                              "complete": complete, "stable_output": stable,
                              "median_seconds": statistics.median(r["wall_seconds"] for r in good) if eligible else "",
                              "min_seconds": min(r["wall_seconds"] for r in good) if eligible else "",
                              "max_seconds": max(r["wall_seconds"] for r in good) if eligible else "",
                              "pair_count": good[0].get("pair_count", "") if eligible else "",
                              "statuses": ";".join(sorted({r["status"] for r in runs}))})
    return summaries


def write_csv(path, rows, fields=None):
    if not rows and fields is None:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def report(output, trials, manifest, tools, repeats):
    summary = summarize(trials, manifest, tools, repeats)
    write_csv(output / "per_structure.csv", summary)
    lookup = {(r["pdb"], r["tool"]): r for r in summary}
    paired = []
    for baseline in ("rubric", "rubric_noscore"):
        if baseline not in tools:
            continue
        for row in summary:
            ref = lookup[(row["pdb"], baseline)]
            if row["tool"] != baseline and row["median_seconds"] != "" and ref["median_seconds"] != "":
                paired.append({"pdb": row["pdb"], "nucleotides": row["nucleotides"],
                               "tool": row["tool"], "baseline": baseline,
                               "tool_seconds": row["median_seconds"], "rubric_seconds": ref["median_seconds"],
                               "speedup": row["median_seconds"] / ref["median_seconds"]})
    write_csv(output / "speedups.csv", paired,
              ["pdb", "nucleotides", "tool", "baseline", "tool_seconds", "rubric_seconds", "speedup"])
    coverage = []
    for tool in tools:
        rows = [r for r in summary if r["tool"] == tool]
        speedups = [r["speedup"] for r in paired if r["tool"] == tool and r["baseline"] == "rubric"]
        coverage.append({"tool": tool, "structures": len(rows),
                         "complete_stable": sum(r["median_seconds"] != "" for r in rows),
                         "paired_with_rubric": len(speedups),
                         "median_speedup_vs_rubric": statistics.median(speedups) if speedups else ""})
    write_csv(output / "coverage.csv", coverage)
