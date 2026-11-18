"""Parsers and chain-renaming helpers for the annotators' output, shared by the comparison
scripts (run_programs.py, convert_outputs.py) and the runtime benchmark (benchmark_direct.py)."""
import re
import string
from itertools import combinations


RV_LINE = re.compile(
    r"^\s*\S+,\s+(\S+):\s+(-?\d+\w?)\s+(\S)-(\S)\s+(-?\d+\w?)\s+(\S+):\s+(\S)/(\S)\s+(cis|tran)"
)

ALIAS_POOL = string.ascii_uppercase + string.digits + string.ascii_lowercase


class TooManyChains(Exception):
    pass


def _atom_site_tags(lines: list[str]) -> list[str]:
    tags: list[str] = []
    for line in lines:
        if line.startswith("_atom_site."):
            tags.append(line.split(".", 1)[1].strip())
        elif tags and line.startswith(("ATOM", "HETATM")):
            break
    return tags


def _atom_site_scan(lines: list[str]) -> tuple[int, int, list[str], set[str]]:
    tags = _atom_site_tags(lines)
    chains: list[str] = []
    nucleic: set[str] = set()
    auth_i = atom_i = -1
    for line in lines:
        if tags and line.startswith(("ATOM", "HETATM")):
            if auth_i < 0:
                auth_i, atom_i = tags.index("auth_asym_id"), tags.index("label_atom_id")
            tok = line.split()
            if len(tok) != len(tags):
                raise ValueError(f"atom_site row has {len(tok)} tokens, expected {len(tags)}")
            if tok[auth_i] not in chains:
                chains.append(tok[auth_i])
            if tok[atom_i].strip('"') == "C1'":
                nucleic.add(tok[auth_i])
    return auth_i, atom_i, chains, nucleic


def _build_alias(chains: list[str], nucleic: set[str]) -> dict[str, str]:
    used: set[str] = set()
    alias: dict[str, str] = {}
    for group in (sorted(nucleic), [c for c in chains if c not in nucleic]):
        for chain in group:
            if len(chain) == 1 and chain not in used:
                free = chain
            else:
                free = next((ch for ch in ALIAS_POOL if ch not in used), None)
            if free is None:
                if chain in nucleic:
                    raise TooManyChains(f"{len(nucleic)} nucleic chains > pool")
                nucleic_aliases = {alias[c] for c in nucleic}
                free = next(ch for ch in ALIAS_POOL if ch not in nucleic_aliases)
            alias[chain] = free
            used.add(free)
    return alias


def chain_groups(names: list[str], size: int = len(ALIAS_POOL)) -> list[list[str]]:
    if len(names) <= size:
        return [list(names)]
    half = size // 2
    blocks = [names[i:i + half] for i in range(0, len(names), half)]
    groups = [blocks[i] + blocks[j] for i, j in combinations(range(len(blocks)), 2)]
    uncovered = {frozenset(p) for p in combinations(names, 2)}
    for group in groups:
        uncovered -= {frozenset(p) for p in combinations(group, 2)}
    if uncovered:
        raise AssertionError(f"chain grouping left {len(uncovered)} pairs uncovered")
    return groups


def _flip_lw(lw: str | None) -> str | None:
    return lw if lw is None or len(lw) != 3 else lw[0] + lw[2] + lw[1]


def _dedupe_pairs(pairs: list[list]) -> list[list]:
    seen = set()
    out = []
    for c1, n1, c2, n2, lw in pairs:
        if (c2, n2) < (c1, n1):
            key = (c2, n2, c1, n1, _flip_lw(lw))
        else:
            key = (c1, n1, c2, n2, lw)
        if key not in seen:
            seen.add(key)
            out.append([c1, n1, c2, n2, lw])
    return out


def alias_cif_chains(text: str) -> tuple[str, dict[str, str]] | None:
    lines = text.splitlines()
    auth_i, _atom_i, chains, nucleic = _atom_site_scan(lines)
    if auth_i < 0 or all(len(c) == 1 for c in chains):
        return None
    alias = _build_alias(chains, nucleic)
    out = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            tok = line.split()
            if len(tok) > auth_i and tok[auth_i] in alias:
                tok[auth_i] = alias[tok[auth_i]]
                line = " ".join(tok)
        out.append(line)
    reverse = {a: c for c, a in alias.items() if c in nucleic or len(c) == 1}
    return "\n".join(out) + "\n", reverse


def _parse_rnaview_out(text: str, reverse: dict[str, str] | None) -> list[list]:
    pairs, in_bp = [], False
    for line in text.splitlines():
        if line.startswith("BEGIN_base-pair"):
            in_bp = True
            continue
        if line.startswith("END_base-pair"):
            break
        if not in_bp or "stack" in line:
            continue
        m = RV_LINE.match(line)
        if not m:
            continue
        c1, n1, _b1, _b2, n2, c2, e1, e2, ct = m.groups()
        if e1 in "+-" and e2 in "+-":
            lw = "cWW"
        elif e1 in "WHS" and e2 in "WHS":
            lw = ("c" if ct == "cis" else "t") + e1 + e2
        else:
            lw = None
        c1, c2 = c1.rstrip(":"), c2.rstrip(":")
        if reverse is not None:
            c1, c2 = reverse[c1], reverse[c2]
        pairs.append([c1, n1, c2, n2, lw])
    return pairs


def _group_cif(lines: list[str], auth_i: int, alias: dict[str, str]) -> str:
    out = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            tok = line.split()
            if tok[auth_i] not in alias:
                continue
            tok[auth_i] = alias[tok[auth_i]]
            line = " ".join(tok)
        out.append(line)
    return "\n".join(out) + "\n"


BARNABA_LW = {"WC": "cWW", "GU": "cWW"}


def _barnaba_lw(ann: str) -> str | None:
    if ann[:2] in BARNABA_LW:
        return BARNABA_LW[ann[:2]]
    if len(ann) == 3 and ann[0] in "WHS" and ann[1] in "WHS" and ann[2] in "ct":
        return ann[2] + ann[0] + ann[1]
    return None


def _bpnet_num(num: str, ins: str) -> str:
    return num if ins == "?" else f"{num}{ins}"


def _parse_bpnet_out(text: str) -> list[list]:
    pairs = []
    for line in text.splitlines():
        if line.startswith("#"):
            continue
        tok = line.split()
        if len(tok) < 13:
            continue
        n1, c1 = _bpnet_num(tok[1], tok[3]), tok[4]
        for k in range(5, len(tok) - 7, 8):
            n2, c2 = _bpnet_num(tok[k + 1], tok[k + 3]), tok[k + 4]
            ptype = tok[k + 5]
            if ":" not in ptype:
                continue
            e1, e2ct = ptype.split(":", 1)
            e1, e2, ct = e1[:1].upper(), e2ct[:1].upper(), e2ct[-1:].upper()
            if e1 in "WHS" and e2 in "WHS" and ct in "CT":
                lw = ("c" if ct == "C" else "t") + e1 + e2
            else:
                lw = None
            pairs.append([c1, n1, c2, n2, lw])
    return pairs
