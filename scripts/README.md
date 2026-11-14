# Analysis scripts and data

Detailed reference for the scripts in this folder and the archived data they read and write. The quick start is in the [main README](../README.md).

## Rebuild the scoring tables from the feature tables
The scoring tables in `data/reference/scoring_tables/` are derived from the
per-pair / per-residue feature CSVs in `data/reference/high_quality_features/`
(and `low_quality_features/` for the penalty weights). Regenerate them with the
self-contained `data_gen` scripts (no engine needed):

```bash
cd scripts/data_gen
python build_distributions.py               # -> parameter_distributions.json
python build_prosco_distributions.py       # -> prosco_distributions.json
python build_z_tables.py                    # -> z_tables.json
python build_backbone_distributions.py      # -> backbone_prosco_distributions.json
python calculate_penalty_weights.py         # -> penalty_weights.json
```

Run them in this order: `calculate_penalty_weights.py` reads the tables the
others write.

## Regenerate the feature tables from raw structures
The feature CSVs come straight from the RUBRIC C++ engine
(https://github.com/YesselmanLab/RUBRIC). Install the engine version that produced
the data in this repository (`requirements-engine.txt` pins `rubric-rna==0.1.0`),
then re-extract; run everything below from the repo root (`RUBRIC-analysis/`).

```bash
# 1. Install the engine (ready-made packages for macOS, Linux and Windows)
pip install -r requirements-engine.txt

# 2. Re-extract every feature table (the engine reads data/structures/<ID>.cif;
#    an entry that is not there is downloaded into data/structures/engine_cache/)
cd scripts/data_gen
python extract_pair_features.py        # -> pairs.csv, hbonds.csv, torsions.csv
python extract_backbone_torsions.py    # -> torsions_all.csv
```

The commands above regenerate the curated reference set. The penalty weights also
need the non-curated set's pair/hbond features:

```bash
python extract_pair_features.py \
    --reference-set ../../data/reference/reference_sets/reference_set_poor.json \
    --output-dir    ../../data/reference/low_quality_features
```

Then rebuild the tables (the rebuild step above) and the figures (the main README). Both extractors take
`--limit N` for a quick smoke test. The scripts use the installed `rubric-rna`;
`--engine /path/to/rubric` or `RUBRIC_ENGINE` points them at another build instead.

## Validation data (PDB-REDO)
The PDB-REDO comparison behind the validation figure is regenerated the same
way - fetch the re-refined coordinates, then rescore original vs REDO with the
engine:

```bash
cd scripts/validation             # from the repo root
python fetch_pdb_redo.py          # uniqueRNAS.csv -> data/validation/redo_files/
python compare_pdb_redo_cpp.py    # score original (RCSB) vs REDO -> the two JSONs
```

`fetch_pdb_redo.py` downloads `<id>_final.pdb` from pdb-redo.eu for every
candidate in `reference_sets/uniqueRNAS.csv` that has a PDB-REDO entry (the rest
are skipped - the ones that succeed *are* the validation set).
`compare_pdb_redo_cpp.py` then writes `pdb_redo_comparison_all.json` and
`pdb_redo_pairs_comparison.json`. PDB-REDO and RCSB are living databases, so a
re-fetch may add or refresh a few entries versus the archived comparison.

## Other analysis scripts

Run from the repo root. Scripts marked *engine* need the RUBRIC engine (the engine step above);
those marked *structures* read `data/structures/` (see [Structures](#structures-datastructures)).

| script | writes | needs |
|---|---|---|
| `scripts/reference/score_pairs_uniques.py` | `data/reference/scores/pair_scores_uniques_rubric.jsonl` | engine, structures |
| `scripts/data_gen/extract_pucker.py` | `data/reference/high_quality_features/pucker.csv` (Fig. S2D) | engine, structures |
| `scripts/data_gen/extract_figure3_examples.py` | `data/reference/figure3/example_pairs.json` (Fig. 3) | engine, structures |
| `scripts/validation/build_redo_suite_points.py` | suite and pucker columns of `data/reference/resolved_again/redo_per_structure_points.csv` (Fig. 4F-H); run before `build_review_points.py` | engine, structures |
| `scripts/validation/build_review_points.py` | `data/reference/resolved_again/redo_per_structure_points.csv`, `issue_rates_per_unit.csv` (Fig. 4, Fig. S1) | engine, structures |
| `scripts/reference/count_nucleotides.py` | nucleotide counts in `data/reference/metadata/metadata_cache.json` | structures |
| `scripts/validation/build_engine_archive.py` | `data/reference/engine_archive/` (per-pair and per-suite engine output) | engine, structures |
| `scripts/validation/build_resolved_cohort.py` | `data/reference/resolved_again/resolve_strict.csv` | |
| `scripts/validation/build_resolved_points.py` | `data/reference/resolved_again/resolve_per_molecule_points.csv` (Fig. 4A-D) | |
| `scripts/validation/build_quality_axes.py` | `data/reference/quality_axes/bang_vs_quality.csv` (Fig. S2B, C) | |
| `scripts/validation/time_reversal_resolved.py` | `resolve_pairs_matched.csv`, `data/validation/time_reversal_resolved.json` | engine |
| `scripts/validation/time_reversal_controls.py` | `data/validation/time_reversal_controls.json` | |
| `scripts/validation/extract_rubric_candidates.py` | `data/comparison/rubric_candidates.json` | engine, structures |
| `scripts/validation/extract_residue_names.py` | `data/comparison/residue_names.json` | structures |
| `scripts/validation/adjudicate_consensus.py` | `data/validation/annotator_consensus.json` | |
| `scripts/validation/consensus_reported_conflicts.py` | `data/validation/annotator_consensus_reported.json` | |
| `scripts/validation/consensus_shared_conflicts.py` | `data/validation/annotator_consensus_shared.json` | |
| `scripts/validation/template_fit_class_conflicts.py` | `data/validation/template_fit_class_conflicts.csv` (Fig. S3D) | structures |
| `scripts/validation/annotator_agreement.py` | `data/validation/annotator_agreement.json` | |
| `scripts/validation/referee_coverage.py` | `data/validation/referee_coverage_current.json` | |
| `scripts/validation/rejection_referee_support.py` | `data/validation/rejection_referee_support.json` | |

The tool-comparison scripts run in the order listed: `adjudicate_consensus.py`
first, then the scripts that read `annotator_consensus.json`.

## Runtime benchmark

The timing of RUBRIC and the six other annotators is documented in
[`benchmark/README.md`](benchmark/README.md); it reads
`data/structures/` and writes `data/benchmark/runtime_direct_all/`.

## Structures (`data/structures/`)

The mmCIF files are not in the archive (12 GB); they are downloaded straight from
RCSB. Every script reads structures from here (the engine's downloads go to
`data/structures/engine_cache/`, see the engine step above), so nothing is read from or written
to anywhere outside the repository. Fetch them with:

```bash
python scripts/fetch_structures.py            # the 2,697 non-redundant RNAs (~12 GB)
python scripts/fetch_structures.py --all      # every RNA entry in the PDB (~10k)
```

The non-redundant set is the `unique_pdb_id` column of
`data/reference/reference_sets/uniqueRNAS.csv`; `--all` uses
`all_pdb_rnas.csv`, which `scripts/fetch_all_pdb_rnas.py` regenerates from the RCSB
search API. Existing files are skipped, so an interrupted download can be repeated.
RCSB occasionally re-issues an entry, so a new download can differ slightly from
the copy used in the paper; the runtime benchmark records each file's checksum.

## What's in `data/reference/`

| subdir | contents |
|---|---|
| `reference_sets/` | the curated reference set, the non-curated set, the non-redundant list (`uniqueRNAS.csv`) and the RCSB RNA universe |
| `high_quality_features/` | per-pair / per-hbond / per-residue features (curated set) |
| `low_quality_features/` | same, for the non-curated set (penalty-weight derivation) |
| `scoring_tables/` | ProSco / Z′ / backbone / penalty-weight tables |
| `scores/`, `metadata/` | per-structure RUBRIC scores and PDB metadata; per-pair scores (`pair_scores_uniques_rubric.jsonl`, one row per selected pair, each non-redundant entry scored once) are regenerated by `scripts/reference/score_pairs_uniques.py`; nucleotide counts in `metadata_cache.json` are recomputed from the local structure files by `scripts/reference/count_nucleotides.py` (counts every residue of an RNA/DNA polymer entity, so modified, locked and L-nucleotides are included) |
| `resolved_again/` | the better-resolution cohort (`resolve_strict.csv`, built by `scripts/validation/build_resolved_cohort.py`) and the per-molecule / per-structure tables behind Figure 4; see its `PROVENANCE.txt` |
| `engine_archive/` | engine output the better-resolution cohort builder reads: per-residue and per-pair rows for the non-redundant set and for the extra members of each BGSU equivalence class (read through `scripts/validation/_archive.py`) |

## PyMOL sessions and images (`data/pymol/`)

| subdir | contents |
|---|---|
| `sessions/` | every PyMOL session (`.pse`) behind the structure panels, each saved in its final camera view; `README.txt` maps sessions to panels |
| `renders/` | the ray-traced images the figure notebooks read (Figures 3 and 4, Supplementary Figures S8, S9 and S12), plus their label anchors (`<key>_meta.json` for Figures S8 and S9, `7D_labels.json` for Figure S12D) |

The notebooks read only `renders/`, never the sessions or the structure files. The
Figure S8 and S9 pairs are rendered from `data/structures` by
`scripts/pymol_figures/render_pair_examples.py`; the other images are rendered from
the sessions by a local script that is not part of the repository.
