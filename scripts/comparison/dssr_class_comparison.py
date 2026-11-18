"""Class-by-class comparison of RUBRIC with DSSR (data for Fig. 2D/E and the conflicts of Fig. S3C/D).

For one LW class and one structure, DSSR's pairs of that class are matched with RUBRIC's
selected pairs of that class. Pairs are restricted to standard A/C/G/U nucleotides; pairs
with a DNA residue and pairs of sequence neighbours are left out. Each unmatched DSSR pair
gets the reason RUBRIC did not report it in this class (from RUBRIC's selection
disposition, or the class RUBRIC assigned instead), and each unmatched RUBRIC pair gets the
class DSSR assigned, if any.
"""
import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

NONCWW_CLASSES = ("tWW", "cWH", "tWH", "cWS", "tWS", "cHH", "tHH", "cHS", "tHS", "cSS", "tSS")
STANDARD_LW = frozenset({"cWW", "cWH", "cWS", "cHH", "cHS", "cSS",
                         "tWW", "tWH", "tWS", "tHH", "tHS", "tSS"})
STANDARD_BASES = {"A", "C", "G", "U"}

BASE_HBOND_ATOMS = frozenset({"N1", "N2", "N3", "N4", "N6", "N7", "N9", "O2", "O4", "O6"})
SUGAR_ATOMS = frozenset({"C1'", "C2'", "C3'", "C4'", "C5'", "O2'", "O3'", "O4'", "O5'"})
PHOSPHATE_ATOMS = frozenset({"P", "OP1", "OP2", "O1P", "O2P", "O5'"})
PAIR_CATEGORY_PRIORITY = ["base-base", "base-sugar", "base-phosphate",
                          "sugar-sugar", "sugar-phosphate", "phosphate-phosphate"]
_HBOND_PATTERN = re.compile(r"([A-Z0-9']+)(?:\([^)]*\))?[-*]([A-Z0-9']+)(?:\([^)]*\))?\[")

# Reason given to an unmatched DSSR pair from RUBRIC's selection disposition of the same
# residue pair. A pair RUBRIC found but did not select keeps 'not_selected'; it becomes
# 'python_assigns:<class>' when RUBRIC classified it in another class.
DISPOSITION_REASON = {
    "selected": "not_selected",
    "competition": "not_selected",
    "ineligible": "no_strong_hbond",
    "stacking": "stacking",
    "adjacent": "not_selected",
    "low_score": "low_score",
    "invalid": "failed_validation",
}
IDENTIFIED_REASONS = {"not_selected", "no_strong_hbond", "stacking", "low_score"}


@dataclass
class DSSRPair:
    res_id1: str
    res_id2: str
    lw_class: str
    saenger: str
    sequence: str
    hbonds_desc: str = ""
    interbase_angle: float = 0.0
    n1n9_distance: float = 0.0
    pair_category: str | None = None
    hbond_categories: list[str] | None = None


@dataclass
class FinderPair:
    res_id1: str
    res_id2: str
    sequence: str
    lw_class: str = "unknown"
    score: float = 0.0
    pair_category: str | None = None
    hbond_categories: list[str] | None = None


@dataclass
class PairDiff:
    index: int = 0
    pdb_id: str = ""
    res_id1: str = ""
    res_id2: str = ""
    sequence: str = ""
    source: str = ""
    reason: str = ""
    pair_score: float | None = None
    dssr_hbonds: str | None = None
    python_hbonds: str | None = None
    template_rmsd: float | None = None
    template_seq: str | None = None
    template_atoms: int | None = None
    dorg: float | None = None
    d_v: float | None = None
    plane_angle: float | None = None
    dNN: float | None = None
    res1_candidates: list | None = None
    res2_candidates: list | None = None
    selected_pairs: list | None = None
    is_helix: bool | None = None
    dssr_lw_class: str | None = None
    python_lw_class: str | None = None
    python_pair_category: str | None = None
    python_lw_class_confidence: list[float] | None = None
    python_lw_class_display: str | None = None
    python_is_ambiguous: bool | None = None
    verdict: str = ""
    verdict_reason: str = ""
    verdict_action: str = ""


@dataclass
class PDBResult:
    pdb_id: str
    dssr_count: int
    python_count: int
    matched: int
    python_only: int
    dssr_only: int
    diffs: list[PairDiff]
    error: str | None = None
    matched_by_seq: dict[str, int] = field(default_factory=dict)


def canonical_lw(lw_class):
    """LW class with its two edges in W < H < S order (cHW -> cWH)."""
    if not lw_class or len(lw_class) < 3:
        return lw_class
    order = {"W": 0, "H": 1, "S": 2}
    if order.get(lw_class[1], 0) > order.get(lw_class[2], 0):
        return lw_class[0] + lw_class[2] + lw_class[1]
    return lw_class


def dssr_to_res_id(dssr_nt):
    """DSSR nucleotide ID to the engine's residue ID ('A.G1^A' -> 'A-G-1A')."""
    if not dssr_nt or "." not in dssr_nt:
        return None
    match = re.match(r"^([A-Za-z0-9]+)\.([A-Za-z0-9]*[A-Za-z][A-Za-z0-9]*?)"
                     r"(-?\d+)([A-Za-z]?)(?:\^([A-Za-z]))?$", dssr_nt)
    if not match:
        return None
    chain, base, num, ins_code, alt_conf = match.groups()
    return f"{chain}-{base}-{num}{ins_code or ''}{alt_conf or ''}"


def _split_res_id(res_id):
    parts = res_id.split("-")
    if len(parts) < 3:
        return None
    return parts[0], parts[1], "-".join(parts[2:])


def are_adjacent(res_id1, res_id2):
    """Same chain and residue numbers differing by at most 1."""
    p1, p2 = _split_res_id(res_id1), _split_res_id(res_id2)
    if p1 is None or p2 is None or p1[0] != p2[0]:
        return False
    try:
        n1 = int("".join(c for c in p1[2] if c.isdigit() or c == "-"))
        n2 = int("".join(c for c in p2[2] if c.isdigit() or c == "-"))
    except ValueError:
        return False
    return abs(n1 - n2) <= 1


def are_adjacent_residues(res_id1, res_id2):
    """Same chain and leading residue numbers differing by at most 1."""
    def number(res_id):
        parts = res_id.split("-")
        digits = ""
        for c in parts[2]:
            if not (c.isdigit() or c == "-"):
                break
            digits += c
        return parts[0], int(digits) if digits else 0
    try:
        (c1, n1), (c2, n2) = number(res_id1), number(res_id2)
    except (ValueError, IndexError):
        return False
    return c1 == c2 and abs(n1 - n2) <= 1


def has_dna_residue(res_id1, res_id2):
    for rid in (res_id1, res_id2):
        name = rid.split("-")[1] if "-" in rid else ""
        if len(name) == 2 and name[0] == "D" and name[1] in "ACGTU":
            return True
    return False


def is_standard_pair(res_id1, res_id2):
    for rid in (res_id1, res_id2):
        parts = rid.split("-")
        if len(parts) < 2 or parts[1] not in STANDARD_BASES:
            return False
    return True


def _atom_group(atom_name):
    name = atom_name.rstrip("*")
    if name in BASE_HBOND_ATOMS:
        return "base"
    if name in SUGAR_ATOMS:
        return "sugar"
    if name in PHOSPHATE_ATOMS:
        return "phosphate"
    return "unknown"


def dssr_pair_category(hbonds_desc):
    """Dominant H-bond category of a DSSR pair from its H-bond description, if given."""
    cats = ["-".join(sorted([_atom_group(m.group(1)), _atom_group(m.group(2))]))
            for m in _HBOND_PATTERN.finditer(hbonds_desc or "")]
    if not cats:
        return "unknown", []
    for category in PAIR_CATEGORY_PRIORITY:
        if category in cats:
            return category, cats
    return cats[0], cats


def parse_dssr_pair(pair_data, lw_class, with_category=True):
    res_id1 = dssr_to_res_id(pair_data.get("nt1", ""))
    res_id2 = dssr_to_res_id(pair_data.get("nt2", ""))
    if not res_id1 or not res_id2:
        return None
    pair = DSSRPair(res_id1=res_id1, res_id2=res_id2, lw_class=lw_class,
                    saenger=pair_data.get("Saenger", ""),
                    sequence=pair_data.get("bp", "").replace("-", "").replace("+", ""),
                    hbonds_desc=pair_data.get("hbonds_desc", ""),
                    interbase_angle=pair_data.get("interBase_angle", 0.0),
                    n1n9_distance=pair_data.get("N1N9_dist", 0.0))
    if with_category:
        pair.pair_category, pair.hbond_categories = dssr_pair_category(pair.hbonds_desc)
        # A pair in one of the twelve LW classes is base-base by definition.
        if canonical_lw(lw_class) in STANDARD_LW:
            pair.pair_category = "base-base"
    return pair


def load_dssr_class(dssr_json, lw_class):
    """DSSR pairs of one LW class (either edge order), keyed by residue pair."""
    target = canonical_lw(lw_class)
    out = {}
    for pair in json.loads(Path(dssr_json).read_text()).get("pairs", []):
        lw = pair.get("LW", "")
        if lw != lw_class and canonical_lw(lw) != target:
            continue
        parsed = parse_dssr_pair(pair, lw)
        if parsed is not None:
            out[(parsed.res_id1, parsed.res_id2)] = parsed
    return out


def load_dssr_all(dssr_json):
    """All DSSR pairs with an LW label, keyed by residue pair."""
    out = {}
    for pair in json.loads(Path(dssr_json).read_text()).get("pairs", []):
        lw = pair.get("LW", "")
        if not lw:
            continue
        parsed = parse_dssr_pair(pair, lw, with_category=False)
        if parsed is not None:
            out[(parsed.res_id1, parsed.res_id2)] = parsed
    return out


def match_pairs(finder_pairs, dssr_pairs):
    """(matched, finder_only, dssr_only), matching residue pairs in either order."""
    matched, finder_only, used = [], [], set()
    for fp in finder_pairs:
        key, rev = (fp.res_id1, fp.res_id2), (fp.res_id2, fp.res_id1)
        if key in dssr_pairs:
            matched.append((fp, dssr_pairs[key]))
            used.add(key)
        elif rev in dssr_pairs:
            matched.append((fp, dssr_pairs[rev]))
            used.add(rev)
        else:
            finder_only.append(fp)
    dssr_only = [p for k, p in dssr_pairs.items() if k not in used]
    return matched, finder_only, dssr_only


def dssr_only_reason(dssr_pair, candidates, selected_keys, min_score=0.1):
    key, rev = (dssr_pair.res_id1, dssr_pair.res_id2), (dssr_pair.res_id2, dssr_pair.res_id1)
    by_key = {(c.res_id1, c.res_id2): c for c in candidates}
    candidate = by_key.get(key) or by_key.get(rev)
    if candidate is None:
        return "not_candidate"
    if candidate.disposition is not None:
        return DISPOSITION_REASON.get(candidate.disposition, "unknown")
    if candidate.quality_score < min_score:
        return "low_score"
    if key not in selected_keys and rev not in selected_keys:
        return "not_selected"
    return "unknown"


def rubric_only_reason(finder_pair, all_dssr):
    dssr_pair = (all_dssr.get((finder_pair.res_id1, finder_pair.res_id2))
                 or all_dssr.get((finder_pair.res_id2, finder_pair.res_id1)))
    if dssr_pair is None:
        return "novel", None
    return f"dssr_assigns:{dssr_pair.lw_class}", dssr_pair.lw_class


def compare_structure(pdb_id, lw_class, dssr_json, selected, candidates, all_dssr):
    """Compare RUBRIC with DSSR for one LW class in one structure.

    selected: RUBRIC's selected pairs as dicts {r1, r2, lw, seq, cat}.
    candidates: RUBRIC's classified candidates (res_id1, res_id2, lw_class, quality_score,
    disposition).
    """
    res = PDBResult(pdb_id=pdb_id, dssr_count=0, python_count=0, matched=0,
                    python_only=0, dssr_only=0, diffs=[])
    try:
        dssr_all_class = load_dssr_class(dssr_json, lw_class)
    except (OSError, ValueError):
        return res
    dssr_pairs = {k: p for k, p in dssr_all_class.items()
                  if p.pair_category == "base-base"
                  and not are_adjacent(p.res_id1, p.res_id2)
                  and not has_dna_residue(p.res_id1, p.res_id2)
                  and is_standard_pair(p.res_id1, p.res_id2)}
    if not dssr_pairs:
        return res

    target = canonical_lw(lw_class)
    rubric_pairs = [FinderPair(res_id1=s["r1"], res_id2=s["r2"], sequence=s["seq"],
                               lw_class=s["lw"] or lw_class)
                    for s in selected
                    if (s["lw"] == lw_class or canonical_lw(s["lw"]) == target)
                    and s["cat"] == "base-base"
                    and not has_dna_residue(s["r1"], s["r2"])
                    and is_standard_pair(s["r1"], s["r2"])]
    selected_lw = {(s["r1"], s["r2"]): s["lw"] for s in selected}
    cand_by_key = {}
    for c in candidates:
        cand_by_key[(c.res_id1, c.res_id2)] = c
        cand_by_key[(c.res_id2, c.res_id1)] = c

    matched, rubric_only, dssr_only = match_pairs(rubric_pairs, dssr_pairs)
    res.dssr_count, res.python_count = len(dssr_pairs), len(rubric_pairs)
    res.matched, res.python_only, res.dssr_only = len(matched), len(rubric_only), len(dssr_only)
    res.matched_by_seq = dict(Counter(dp.sequence for _, dp in matched))

    selected_keys = {(s["r1"], s["r2"]) for s in selected}
    for dp in dssr_only:
        if are_adjacent_residues(dp.res_id1, dp.res_id2):
            continue
        reason = dssr_only_reason(dp, candidates, selected_keys)
        if reason in ("not_selected", "unknown", "no_strong_hbond", "stacking"):
            cand = cand_by_key.get((dp.res_id1, dp.res_id2))
            own_lw = (selected_lw.get((dp.res_id1, dp.res_id2))
                      or selected_lw.get((dp.res_id2, dp.res_id1)))
            if own_lw is None and cand is not None:
                own_lw = cand.lw_class
            score = cand.quality_score if cand else 0.0
            if own_lw and score >= 0.5 and own_lw != lw_class and canonical_lw(own_lw) != target:
                reason = f"python_assigns:{own_lw}"
        res.diffs.append(PairDiff(pdb_id=pdb_id, res_id1=dp.res_id1, res_id2=dp.res_id2,
                                  sequence=dp.sequence, source="dssr_only",
                                  reason=reason, is_helix=False))

    for fp in rubric_only:
        reason, dssr_lw = rubric_only_reason(fp, all_dssr)
        if dssr_lw and "." in dssr_lw:
            res.python_only -= 1
        res.diffs.append(PairDiff(pdb_id=pdb_id, res_id1=fp.res_id1, res_id2=fp.res_id2,
                                  sequence=fp.sequence, source="python_only",
                                  reason=reason, dssr_lw_class=dssr_lw))
    return res


def summarize(results, lw_class):
    """Totals, recall and precision of one LW class over all structures."""
    total_dssr = sum(r.dssr_count for r in results)
    total_rubric = sum(r.python_count for r in results)
    total_matched = sum(r.matched for r in results)
    matched_by_sequence = defaultdict(int)
    for r in results:
        for seq, count in r.matched_by_seq.items():
            matched_by_sequence[seq] += count
    dssr_only = [d for r in results for d in r.diffs if d.source == "dssr_only"]
    rubric_only = [d for r in results for d in r.diffs if d.source == "python_only"]

    by_reason_dssr, by_sequence_dssr, by_reason_helix = Counter(), defaultdict(int), Counter()
    helix = isolated = 0
    for d in dssr_only:
        by_reason_dssr[d.reason] += 1
        by_sequence_dssr[d.sequence] += 1
        by_reason_helix[(d.reason, bool(d.is_helix))] += 1
        helix += bool(d.is_helix)
        isolated += not d.is_helix
    by_reason_rubric, by_dssr_class = Counter(), Counter()
    for d in rubric_only:
        by_reason_rubric[d.reason] += 1
        if d.dssr_lw_class:
            by_dssr_class[d.dssr_lw_class] += 1

    target = canonical_lw(lw_class)
    classified_correct = total_matched
    for d in dssr_only:
        if d.reason in IDENTIFIED_REASONS:
            classified_correct += 1
        elif d.reason.startswith("python_assigns:"):
            assigned = d.reason.split(":", 1)[1]
            if "|" in assigned and target in [canonical_lw(p) for p in assigned.split("|")]:
                classified_correct += 1
    misclassified = novel = 0
    for d in rubric_only:
        if d.reason.startswith("dssr_assigns:"):
            misclassified += canonical_lw(d.reason.split(":", 1)[1]) != target
        elif d.reason == "novel":
            novel += 1

    ratio = lambda a, b: a / b if b > 0 else 0
    recall, precision = ratio(total_matched, total_dssr), ratio(total_matched, total_rubric)
    f1 = ratio(2 * recall * precision, recall + precision)
    pdb_recalls = [r.matched / r.dssr_count for r in results if r.dssr_count > 0]
    pdb_precisions = [r.matched / r.python_count for r in results if r.python_count > 0]
    return {
        "lw_class": lw_class,
        "total_pdbs": len(results),
        "errors": sum(1 for r in results if r.error),
        "total_dssr_pairs": total_dssr,
        "total_python_pairs": total_rubric,
        "matched": total_matched,
        "python_only": sum(r.python_only for r in results),
        "dssr_only": sum(r.dssr_only for r in results),
        "recall": round(recall, 4),
        "precision": round(precision, 4),
        "macro_recall": round(ratio(sum(pdb_recalls), len(pdb_recalls)), 4),
        "macro_precision": round(ratio(sum(pdb_precisions), len(pdb_precisions)), 4),
        "n_pdbs_with_dssr": len(pdb_recalls),
        "n_pdbs_with_python": len(pdb_precisions),
        "f1_score": round(f1, 4),
        "dssr_only_helix": helix,
        "dssr_only_isolated": isolated,
        "dssr_only_by_reason": dict(by_reason_dssr),
        "dssr_only_by_reason_context": {
            reason: {"total": n, "helix": by_reason_helix.get((reason, True), 0),
                     "isolated": by_reason_helix.get((reason, False), 0)}
            for reason, n in by_reason_dssr.items()},
        "dssr_only_by_sequence": dict(by_sequence_dssr),
        "matched_by_sequence": dict(matched_by_sequence),
        "python_only_by_reason": dict(by_reason_rubric),
        "python_only_dssr_classes": dict(by_dssr_class),
        "classified_correct": classified_correct,
        "classified_recall": round(ratio(classified_correct, total_dssr), 4),
        "classified_misclassified": misclassified,
        "classified_precision": round(ratio(classified_correct, classified_correct + misclassified), 4),
        "novel_pairs": novel,
    }
