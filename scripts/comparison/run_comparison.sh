#!/usr/bin/env bash
# Regenerates the comparison with the six annotation programs, from program runs to the
# files behind Fig. 2, Fig. S3, Fig. S9 and Table S5. Run from the repository root.
# Set SKIP_RUNS=1 to start from the kept program outputs in data/comparison/program_runs
# (the archived outputs) instead of running the programs.
set -euo pipefail
cd "$(dirname "$0")/../.."
PY="${PYTHON:-python}"
TOOLS_PY="${TOOLS_PYTHON:-.runtime-tools/venv/bin/python}"

if [ "${SKIP_RUNS:-0}" != "1" ]; then
    "$TOOLS_PY" scripts/comparison/run_programs.py
    "$PY" scripts/comparison/keep_rubric_output.py
fi
"$PY" scripts/comparison/convert_outputs.py
"$PY" scripts/comparison/extract_rubric.py
"$PY" scripts/comparison/build_comparison_set.py
"$PY" scripts/comparison/dssr_json.py
"$PY" scripts/comparison/build_per_class.py
"$PY" scripts/comparison/annotate_dssr_candidates.py

for s in adjudicate_consensus consensus_shared_conflicts template_fit_class_conflicts \
         tool_support_profile; do
    "$PY" "scripts/comparison/$s.py"
done
