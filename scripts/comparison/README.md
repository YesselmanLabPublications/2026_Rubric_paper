# Comparison with six annotation programs

These scripts regenerate the data behind Fig. 2B-F, Fig. S3, Fig. S9 and Table S5. RUBRIC is compared with DSSR, FR3D, RNAView, barnaba, RNApolis Annotator and BPNet on the non-redundant set (`data/reference/reference_sets/uniqueRNAS.csv`, 2,696 entries).

## Requirements
- **The six programs.** They are installed with `scripts/benchmark/install_tools.py` into `.runtime-tools/`, with the versions and commands of Table S4 (see [`../benchmark/README.md`](../benchmark/README.md)). DSSR needs a licence from its authors and is not redistributed.
- **RUBRIC 1.0**, installed with `pip install -r requirements-engine.txt`.
- **The structure files** in `data/structures/`, which come with `rubric_structures.zip` (`bash scripts/fetch_data.sh structures`).

## Reproducing without running the programs
The output of every program run, and RUBRIC's own output, come with `rubric_program_runs.zip`
(`bash scripts/fetch_data.sh program_runs`), in `data/comparison/program_runs/`. Start from those outputs:

```bash
SKIP_RUNS=1 bash scripts/comparison/run_comparison.sh
```

This also reproduces the DSSR-based results without a DSSR licence.

## Full run
```bash
bash scripts/comparison/run_comparison.sh
```

Then re-run `notebooks/figures/figure2.ipynb`, `notebooks/supplemental_figures/figS3.ipynb` and `notebooks/supplemental_figures/figS9.ipynb`.

| Step | Script | Writes |
|---|---|---|
| 1 | `run_programs.py`, `keep_rubric_output.py` | `data/comparison/program_runs/raw/`, `status.jsonl` (native output of every run; RUBRIC's JSON output in `raw/rubric/`) |
| 2 | `convert_outputs.py` | `program_runs/pairs/<program>_pairs.json` |
| 3 | `extract_rubric.py` | `program_runs/rubric/` (RUBRIC pairs, candidates, residue identities) |
| 4 | `build_comparison_set.py` | `data/comparison/<program>_pairs.json`, `rubric_*.json`, `residue_names.json`, `fr3d/`, `comparison_entries.json` |
| 5 | `dssr_json.py` | `program_runs/dssr_json/` |
| 6 | `build_per_class.py` | `data/comparison/<class>.json` |
| 7 | `annotate_dssr_candidates.py` | adds RUBRIC's class to each DSSR pair RUBRIC did not report |
| 8 | `adjudicate_consensus.py`, `consensus_shared_conflicts.py`, `template_fit_class_conflicts.py`, `tool_support_profile.py` | `data/validation/annotator_consensus*.json`, `template_fit_class_conflicts.csv`, `tool_support_profile.json` (Fig. S3) |

`comparison_data.py` loads the comparison files for the Fig. S3 notebook and `template_fit_class_conflicts.py`.

`run_programs.py` uses the Python environment of the programs (`.runtime-tools/venv`). The other steps use the analysis environment.

## Rules that affect the results
- **Commands.** Each program runs with the command of Table S4 on the deposited mmCIF file.
- **RNAView chain IDs.** RNAView reads only the first character of a chain ID. For entries with longer chain IDs it runs on a copy with the chains renamed to single characters, and the names are restored in its output. Entries with more nucleic-acid chains than single characters run in groups of chains that together cover every pair of chains.
- **Which outputs are used.** An output counts if the program exited normally, or if it exited with an error after writing a complete base-pair list. For RNAView, a complete list ends with `END_base-pair`. For BPNet, it lists every residue that its header counts.
- **Rows left out:**
  - rows that pair a residue with itself (RNAView and barnaba print residue numbers without insertion codes);
  - FR3D rows that involve a symmetry-related copy, or an alternate location other than A.
- **Agreement statistics** (Fig. 2B, C, F; Fig. S3; Fig. S9; Table S5) use the entries that RUBRIC and all six programs processed (`comparison_entries.json`).
- **Recall and precision by class** (Fig. 2D/E) use every entry with DSSR output. They are restricted to pairs of standard A, C, G and U nucleotides, without DNA residues and without sequence neighbours. A DSSR pair in one of the twelve Leontis-Westhof classes counts as a base-base pair.
