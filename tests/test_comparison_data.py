"""Hand-counted cases for the comparison with the annotation programs."""

import json

import pytest

from comparison_data import Comparison, flip_label, ingest, residue

SELECTED = [["A", "1", "A", "4", "cWW"], ["A", "40", "A", "30", "cSH"]]
CANDIDATES = [
    [*SELECTED[0], "selected"],
    [*SELECTED[1], "selected"],
    ["A", "10", "A", "20", "tWH", "competition"],
    ["A", "50", "A", "60", "tWW", "ineligible"],
    ["A", "70", "A", "80", "cWW|tWW", "stacking"],
    ["A", "90", "A", "100", None, "ineligible"],
    ["A", "110", "A", "120", "cWW", "invalid"],
]
DSSR = [
    ["A", "1", "A", "4", "cWW"],
    ["A", "10", "A", "20", "tWH"],
    ["A", "30", "A", "40", "cWH"],
    ["A", "70", "A", "80", "tWW"],
]


@pytest.fixture
def comparison(tmp_path):
    folder = tmp_path / "comparison"
    folder.mkdir()
    files = {
        "rubric_pairs.json": {p: {"status": "ok", "pairs": SELECTED} for p in ["P1", "OLD"]},
        "rubric_candidates.json": {
            "P1": {"status": "ok", "pairs": CANDIDATES},
            "OLD": {"status": "ok", "pairs": CANDIDATES},
            "_meta": {"label_order": "input_residue_order", "fallback_to_output": ["OLD"]},
        },
        "dssr_pairs.json": {"P1": {"status": "ok", "pairs": DSSR}},
        "barnaba_pairs.json": {"P1": {"status": "ok", "pairs": []}},
        "residue_names.json": {"P1": {"A|1|": "G", "A|4|": "C"}},
    }
    for name, value in files.items():
        (folder / name).write_text(json.dumps(value))
    return Comparison(tmp_path)


def test_identification_does_not_depend_on_comparator(comparison):
    original = {k: set(v) for k, v in comparison.identified["P1"].items()}
    comparison.detection_counts("rubric", "dssr")
    comparison.detection_counts("rubric", "barnaba")
    assert comparison.identified["P1"] == original
    assert len(original) == 5
    assert comparison.metrics("identified")["tWW"]["total_python"] == 1
    assert comparison.metrics("identified")["tWW"]["prec_cls"] == 0


def test_credit_requires_the_comparator_to_report_the_pair(comparison):
    credited = comparison.credited("dssr")["P1"]
    assert (residue("A", "10"), residue("A", "20")) in credited
    assert (residue("A", "50"), residue("A", "60")) not in credited
    assert len(credited) == 4
    assert comparison.credited("barnaba")["P1"] == comparison.selected["P1"]
    assert comparison.metrics()["tWW"]["total_python"] == 0
    assert comparison.metrics()["tWH"]["matched"] == 1


def test_credit_keeps_rubric_class_when_the_comparator_disagrees(comparison):
    key = (residue("A", "10"), residue("A", "20"))
    comparison.tool("dssr")["P1"][key] = {"cWH"}
    assert comparison.credited("dssr")["P1"][key] == {"tWH"}
    assert comparison.metrics()["tWH"]["total_python"] == 1
    assert comparison.metrics()["tWH"]["matched"] == 0


def test_selected_pairs_count_without_comparator_support(comparison):
    _, shared, union, _ = comparison.detection_counts("rubric", "barnaba")
    assert shared.sum() == 0
    assert union.sum() == 2


def test_detection_is_symmetric_and_class_counts_match(comparison):
    forward = comparison.detection_counts("rubric", "dssr")
    reverse = comparison.detection_counts("dssr", "rubric")
    assert forward[0] == reverse[0]
    assert (forward[1] == reverse[1]).all()
    assert (forward[2] == reverse[2]).all()
    assert forward[1].sum() == 3
    assert forward[2].sum() == 3
    assert comparison.detection_counts(stage="identified")[2].sum() == 4
    credited = comparison.class_counts()
    identified = comparison.class_counts(stage="identified")
    assert (credited[1] == identified[1]).all()
    assert (credited[2] == identified[2]).all()


def test_same_reference_denominator_at_every_stage(comparison):
    effect = comparison.selection_effect()
    for kind in ["canonical", "non-canonical"]:
        reference = effect.loc[("selected", kind), "reference"]
        assert effect.loc[("identified", kind), "reference"] == reference
        assert effect.loc[("credited", kind), "reference"] == reference
    assert effect.loc[("identified", "non-canonical"), "precision"] == pytest.approx(1 / 3)
    assert effect.loc[("identified", "non-canonical"), "recall"] == pytest.approx(1 / 2)
    assert effect.loc[("selected", "non-canonical"), "recall"] == 0
    assert effect.loc[("credited", "non-canonical"), "precision"] == pytest.approx(1 / 2)
    assert effect.loc[("credited", "non-canonical"), "recall"] == pytest.approx(1 / 2)
    assert comparison.dssr_stage_counts[1] == 1


def test_labels_follow_residue_order_and_insertion_codes():
    a, b = residue("AA", "-2"), residue("AA", "10A")
    assert ingest([(b, a, "cHW")]) == {(a, b): {"cWH"}}
    assert flip_label("cWH|tHS") == "cHW|tSH"
    assert residue("AA", "10", "A") == b


def test_class_disagreement_is_not_a_match(comparison):
    metrics = comparison.metrics()
    assert metrics["cHS"]["matched"] == 0
    assert metrics["cHS"]["total_python"] == 1
    assert metrics["cWH"]["total_dssr"] == 1


def test_mismatched_entries_excluded_from_every_stage(comparison):
    assert comparison.cohort == {"P1"}
    assert comparison.selected.keys() == comparison.identified.keys()
    assert comparison.selected["P1"].keys() <= comparison.identified["P1"].keys()


def test_unclassified_and_invalid_candidates_not_counted(comparison):
    counts = comparison.disposition_counts()
    assert counts.sum() == 5
    assert counts["selected"] == 2
    assert counts["ineligible"] == 1
