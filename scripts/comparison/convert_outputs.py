"""Convert the kept program outputs to base-pair lists.

Writes data/comparison/program_runs/pairs/<program>_pairs.json as
{PDB: {"status": "ok" | "error", "via": ..., "pairs": [[chain1, number1, chain2, number2, lw], ...]}}.
Residue numbers carry their insertion code ("72A"); lw is the program's Leontis-Westhof
label, or None when it gives no single class. FR3D is written as fr3d_local_pairs.json.

An output is used when the program exited normally, or when it exited with an error after
writing its complete base-pair list: RNAView's list ends with END_base-pair, and BPNet's
lists every residue its header counts. Rows are left out when both residues are the same
(RNAView and barnaba print residue numbers without insertion codes, so 129A-129B reads as
129-129), and FR3D rows that involve a symmetry-related copy or an alternate location other
than A, which the other programs do not report.
"""
import argparse
import gzip
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "benchmark"))
from _annotators import (_atom_site_scan, _barnaba_lw, _build_alias, _dedupe_pairs,  # noqa: E402
                         _parse_bpnet_out, _parse_rnaview_out, alias_cif_chains, chain_groups)

RUNS = ROOT / "data" / "comparison" / "program_runs"
LW_CLASS = re.compile(r"^[ct][WHS][WHS]$")


def read_gz(path):
    with gzip.open(path, "rt", errors="replace") as fh:
        return fh.read()


def dssr_residue(nt):
    """'A.G25' -> ('A', '25'); 'A.5MC/40^B' -> ('A', '40B'); 'A.G-5' -> ('A', '-5')."""
    chain, rest = nt.split(".", 1)
    icode = ""
    if "^" in rest:
        rest, icode = rest.split("^", 1)
    if "/" in rest:
        number = rest.rsplit("/", 1)[1]
    else:
        m = re.match(r"^(.*?[A-Za-z_])(-?\d+)$", rest) or re.match(r"^()(-?\d+)$", rest)
        if not m:
            raise ValueError(f"unparsed DSSR residue {nt!r}")
        number = m.group(2)
    return chain, number + icode


def parse_dssr(text):
    """DSSR's base-pair table; the LW column as written ('...' for no class is kept as '--')."""
    rows, inside, expected = [], False, None
    for line in text.splitlines():
        m = re.match(r"^List of (\d+) base pairs?$", line.strip())
        if m:
            inside, expected = True, int(m.group(1))
            continue
        if not inside:
            continue
        tok = line.split()
        if tok and tok[0] == "nt1":
            continue
        if not tok or not tok[0].isdigit():
            if rows:
                break
            continue
        c1, n1 = dssr_residue(tok[1])
        c2, n2 = dssr_residue(tok[2])
        rows.append([c1, n1, c2, n2, "--" if tok[-2] == "..." else tok[-2]])
    if expected is not None and len(rows) != expected:
        raise ValueError(f"DSSR rows {len(rows)} != {expected}")
    return rows


def fr3d_unit_in_model(unit):
    """False for a unit in a symmetry-related copy or at an alternate location other than A."""
    f = unit.split("|")
    altloc = f[6] if len(f) >= 7 else ""
    symop = f[8] if len(f) >= 9 else ""
    return altloc in ("", "A") and symop in ("", "1_555")


def fr3d_residue(unit):
    f = unit.split("|")
    if len(f) < 5 or f[1] != "1":
        return None
    return f[2], f[4] + (f[7] if len(f) >= 8 else "")


def parse_fr3d(text):
    rows = []
    for line in text.splitlines():
        f = line.split("\t")
        if len(f) < 3 or not (fr3d_unit_in_model(f[0]) and fr3d_unit_in_model(f[2])):
            continue
        a, b = fr3d_residue(f[0]), fr3d_residue(f[2])
        if a and b:
            rows.append([a[0], a[1], b[0], b[1], f[1]])
    return rows


def parse_rnapolis(text):
    rows = []
    for bp in json.loads(text)["base_pairs"]:
        def residue(nt):
            a = nt["auth"]
            icode = (a.get("icode") or "").strip()
            return a["chain"], f"{a['number']}{'' if icode in ('', '?') else icode}"
        c1, n1 = residue(bp["nt1"])
        c2, n2 = residue(bp["nt2"])
        lw = bp.get("lw")
        rows.append([c1, n1, c2, n2, lw if lw and LW_CLASS.match(lw) else None])
    return rows


def mdtraj_chain_auth_ids(cif):
    """Author chain ID of each chain mdtraj creates from an mmCIF file, in mdtraj's order.

    mdtraj uses label_asym_id as the chain column when it has more distinct values than
    auth_asym_id, otherwise auth_asym_id, and starts a new chain whenever either column
    changes between consecutive atom rows of the first model.
    """
    import gemmi
    block = gemmi.cif.read(str(cif)).sole_block()
    auth = list(block.find_values("_atom_site.auth_asym_id"))
    label = list(block.find_values("_atom_site.label_asym_id"))
    model = list(block.find_values("_atom_site.pdbx_PDB_model_num")) or ["1"] * len(auth)
    use_label = len(set(label)) > len(set(auth))
    chains, last = [], None
    for a, lab, m in zip(auth, label, model):
        if m != model[0]:
            continue
        key = (lab, a) if use_label else (a, lab)
        if key != last:
            chains.append(a)
            last = key
    return chains


def parse_barnaba(text, cif):
    """barnaba names residues <base>_<number>_<mdtraj chain index>."""
    chain_ids = mdtraj_chain_auth_ids(cif)
    rows = []
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        r1, r2, ann = line.split()[:3]
        _b1, n1, i1 = r1.rsplit("_", 2)
        _b2, n2, i2 = r2.rsplit("_", 2)
        rows.append([chain_ids[int(i1)], n1, chain_ids[int(i2)], n2, _barnaba_lw(ann)])
    return rows


def parse_rnaview(record, runs, cif):
    """RNAView pairs with the original chain IDs restored for renamed and grouped runs."""
    texts = [read_gz(runs / o) for o in record["outputs"]]
    via = record["via"]
    if via == "direct":
        return _parse_rnaview_out(texts[0], None)
    if via == "renamed":
        _text, reverse = alias_cif_chains(cif.read_text(errors="replace"))
        return _parse_rnaview_out(texts[0], reverse)
    size = int(via.rsplit("_", 1)[1])
    lines = cif.read_text(errors="replace").splitlines()
    _auth_i, _atom_i, _chains, nucleic = _atom_site_scan(lines)
    rows = []
    for text, group in zip(texts, chain_groups(sorted(nucleic), size)):
        alias = _build_alias(group, set(group))
        rows += _parse_rnaview_out(text, {a: c for c, a in alias.items()})
    return _dedupe_pairs(rows)


def bpnet_complete(text):
    n = next((int(l.split(":")[1]) for l in text.splitlines() if "Cleaned number of residues" in l), None)
    last = max((int(l.split()[0]) for l in text.splitlines()
                if l.strip() and not l.startswith("#") and l.split()[0].isdigit()), default=0)
    return n is not None and last == n


def complete(record, runs):
    """Whether the kept output can be used (see module docstring)."""
    if not record["outputs"]:
        return False
    if all(c == 0 for c in record["returncodes"]):
        return True
    texts = [read_gz(runs / o) for o in record["outputs"]]
    if record["program"] == "rnaview":
        return all("END_base-pair" in t for t in texts)
    if record["program"] == "bpnet":
        return bpnet_complete(texts[0])
    return False


def convert(record, runs, structures):
    program, cif = record["program"], structures / f"{record['pdb']}.cif"
    if program == "rnaview":
        rows = parse_rnaview(record, runs, cif)
    else:
        text = read_gz(runs / record["outputs"][0])
        rows = {"dssr": parse_dssr, "fr3d": parse_fr3d, "rnapolis": parse_rnapolis,
                "bpnet": _parse_bpnet_out}.get(program, lambda t: parse_barnaba(t, cif))(text)
    return [r for r in rows if not (r[0] == r[2] and str(r[1]) == str(r[3]))]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--runs", type=Path, default=RUNS)
    ap.add_argument("--structures", type=Path, default=ROOT / "data" / "structures")
    args = ap.parse_args()
    latest = {}
    for line in (args.runs / "status.jsonl").read_text().splitlines():
        r = json.loads(line)
        latest[(r["pdb"], r["program"])] = r
    out = args.runs / "pairs"
    out.mkdir(exist_ok=True)
    for program in sorted({p for _pdb, p in latest}):
        result = {}
        for (pdb, prog), record in sorted(latest.items()):
            if prog != program:
                continue
            if not complete(record, args.runs):
                result[pdb] = {"status": "error", "via": record["via"], "pairs": []}
                continue
            result[pdb] = {"status": "ok", "via": record["via"],
                           "pairs": convert(record, args.runs, args.structures)}
        name = "fr3d_local" if program == "fr3d" else program
        (out / f"{name}_pairs.json").write_text(json.dumps(result))
        ok = sum(v["status"] == "ok" for v in result.values())
        print(f"{name}: {ok} of {len(result)} entries usable, "
              f"{sum(len(v['pairs']) for v in result.values()):,} pairs")


if __name__ == "__main__":
    main()
