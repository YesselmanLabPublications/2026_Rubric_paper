"""Geometry of the RUBRIC-DSSR class conflicts (data for Fig. S3D).

Where RUBRIC and DSSR detect the same pair but assign different LW classes, the five
other programs mostly side with RUBRIC (Fig. S3C). Agreement among programs does not
show which class is right, so this script adds a geometric line of evidence: for every
conflict, it asks whether the pair's geometry is closer to the template of RUBRIC's
class or of DSSR's class.

The templates are the ones RUBRIC uses (data/reference/basepair_idealized, a copy of the
engine's resources/basepair-idealized, built by RUBRIC/scripts/build_templates.py). For
every class and sequence, the template is the medoid of pairs from the curated reference
set that the six other programs (DSSR, FR3D, RNAView, barnaba, RNApolis and BPNet) assign
to that class. The conflicts are the population of Fig. S3C: pairs both programs detect,
with a definite DSSR class, that at least one referee program reports.

Each conflict is superposed (Kabsch, base atoms of both residues) on the template of each
class, with the residues in either order (classes are compared without orientation, as in
Fig. S3C), and the better fit is kept. RMSD differences within 0.1 A count as undecided.

Writes
  data/validation/template_fit_class_conflicts.csv      (read by Fig. S3D)
  figures/source_data/template_fit_class_conflicts.csv  (source data of Fig. S3D)
and prints the share of conflicts closer to each class, with 95% bootstrap confidence
intervals (2,000 resamples of structures) and a two-sided z-test with cluster-robust
variance over structures, overall and without the largest conflict type, by conflict
type, and against the referee verdict.

Caveats
- The templates are built from pairs that the six other programs agree on, so RUBRIC's
  own class assignments play no part in them. DSSR is one of the six, so where a conflict
  pair resembles pairs that DSSR and the others agree on, the test favours DSSR.
- RUBRIC uses the same templates to resolve ambiguous pairs, so for those pairs the
  comparison is not independent of how RUBRIC chose the class.
- Conflicts without a template for both classes are left out; the count is printed.
- The largest conflict type (tSS in RUBRIC, tWS in DSSR) dominates the overall share, so
  the result is also given without it.

Needs the mmCIF files in data/structures (see scripts/fetch_structures.py).

Usage (from the repo root):  python scripts/validation/template_fit_class_conflicts.py
"""
import sys
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "validation"))
import consensus_shared_conflicts as conflicts  # noqa: E402
from _vendor import annotators  # noqa: E402
from figure2_comparison import Comparison  # noqa: E402

data_file = lambda rel: ROOT / "data" / rel  # noqa: E731
SOURCE_DATA = ROOT / "figures" / "source_data"
TIE = 0.1                 # RMSD differences within 0.1 A count as undecided
ATOMS = {"A": ["C1'", "N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"],
         "G": ["C1'", "N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"],
         "C": ["C1'", "N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"],
         "U": ["C1'", "N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"]}


def canon(lw):
    order = {"W": 0, "H": 1, "S": 2}
    return lw[0] + lw[2] + lw[1] if order[lw[1]] > order[lw[2]] else lw


def class_conflicts():
    """The Fig. S3C population, with the votes of the referee programs for each side."""
    panels = {"fr3d": annotators.fr3d_pairs(),
              **{t: annotators.tool_pairs(t) for t in conflicts.VOTING_TOOLS}}
    rows = []
    for r in conflicts.class_records(conflicts.s1a_conflicts()):
        found = {t: {c for c in annotators.lookup(p[r.pdb_id], r.res1, r.res2) if not c.startswith("n")}
                 for t, p in panels.items() if r.pdb_id in p}
        sides = annotators.vote_sides(r.disputed_class, r.rubric_candidates)
        if not found or sides is None:
            continue
        dssr_side, rubric_side = sides
        rows.append({"pdb": r.pdb_id, "key": tuple(sorted([(r.res1.chain, r.res1.number, r.res1.icode),
                                                          (r.res2.chain, r.res2.number, r.res2.icode)])),
                     "dssr": r.disputed_class, "rubric": [c for c in r.rubric_candidates if c not in dssr_side],
                     "votes_rubric": sum(bool(c & rubric_side) for c in found.values()),
                     "votes_dssr": sum(bool(c & dssr_side) for c in found.values())})
    return pd.DataFrame(rows)


def find_residue(model, rid):
    chain = model.find_chain(rid[0])
    return next((r for r in chain or [] if r.seqid.num == rid[1] and r.seqid.icode.strip() == rid[2]), None)


def coords(res, base):
    pts = []
    for name in ATOMS[base]:
        atom = res.find_atom(name, "*")
        if atom is None:
            return None
        pts.append(atom.pos.tolist())
    return np.array(pts)


def tier_of(path):
    for line in path.read_text().splitlines():
        if line.startswith("REMARK   2 Tier"):
            return int(line.split()[3].rstrip(":"))
    return None


def load_templates():
    """One template per class and sequence; symmetric classes are stored in both orders."""
    templates, tier = {}, {}
    for path in sorted(data_file("reference/basepair_idealized").glob("*/*.pdb")):
        lw, seq = path.parent.name, path.stem
        if len(seq) != 2 or not set(seq) <= set("ACGU"):
            continue
        residues = [r for ch in gemmi.read_structure(str(path))[0] for r in ch]
        if len(residues) != 2:
            continue
        xa, xb = coords(residues[0], seq[0]), coords(residues[1], seq[1])
        if xa is None or xb is None:
            continue
        templates[(lw, seq)] = np.vstack([xa, xb])
        tier[(lw, seq)] = tier_of(path)
    return templates, tier


def superpose(p, q):
    """Rotate/translate q onto p; return transformed q and RMSD."""
    pc, qc = p.mean(0), q.mean(0)
    u, s, vt = np.linalg.svd((q - qc).T @ (p - pc))
    d = np.sign(np.linalg.det(u @ vt))
    rot = u @ np.diag([1, 1, d]) @ vt
    moved = (q - qc) @ rot + pc
    return moved, float(np.sqrt(((moved - p) ** 2).sum(1).mean()))


def fit(templates, lw, residues, bases):
    """RMSD to the template of class lw, residues in either order; the better fit."""
    best = None
    for (x, bx), (y, by) in (((residues[0], bases[0]), (residues[1], bases[1])),
                              ((residues[1], bases[1]), (residues[0], bases[0]))):
        tmpl = templates.get((canon(lw), bx + by))
        if tmpl is None:
            continue
        cx, cy = coords(x, bx), coords(y, by)
        if cx is None or cy is None:
            continue
        value = superpose(tmpl, np.vstack([cx, cy]))[1]
        best = value if best is None else min(best, value)
    return best


def summarize(df, reps=2000, seed=7):
    per = df.groupby("pdb").closer_to.agg(n="size", a=lambda s: (s == "RUBRIC").sum(),
                                          b=lambda s: (s == "DSSR").sum())
    n, a, b = (per[c].to_numpy(float) for c in ("n", "a", "b"))
    w = np.random.default_rng(seed).multinomial(len(n), np.full(len(n), 1 / len(n)), size=reps).astype(float)
    diff = (a - b).sum() / n.sum()
    u = (a - b - diff * n) / n.sum()
    p = 2 * norm.sf(abs(diff / np.sqrt(len(u) / (len(u) - 1) * (u ** 2).sum())))
    ci = lambda v: np.percentile(100 * (w @ v) / (w @ n), [2.5, 97.5]).round(1)  # noqa: E731
    return {"pairs": int(n.sum()), "structures": len(n),
            "closer to RUBRIC (%)": round(100 * a.sum() / n.sum(), 1), "95% CI RUBRIC": ci(a),
            "closer to DSSR (%)": round(100 * b.sum() / n.sum(), 1), "95% CI DSSR": ci(b),
            "within 0.1 A (%)": round(100 * (n - a - b).sum() / n.sum(), 1),
            "p (cluster-robust)": "< 0.0001" if p < 1e-4 else f"{p:.2g}"}


def main():
    conflict_table = class_conflicts()
    print(f"{len(conflict_table):,} class conflicts in {conflict_table.pdb.nunique():,} structures")

    names = Comparison(ROOT / "data").residue_names
    templates, tier = load_templates()
    print(f"{len(templates)} templates (class and sequence, both orders for symmetric classes); by tier",
          dict(sorted(pd.Series(tier).value_counts().items())))

    conflict_res = {}
    for pdb, group in conflict_table.groupby("pdb"):
        model = gemmi.read_structure(str(data_file(f"structures/{pdb}.cif")))[0]
        for r in group.itertuples():
            residues = [find_residue(model, k) for k in r.key]
            conflict_res[r.Index] = residues if None not in residues else None
    print(f"coordinates for {sum(v is not None for v in conflict_res.values()):,} "
          f"of {len(conflict_table):,} conflicts")

    records = []
    for r in conflict_table.itertuples():
        residues = conflict_res.get(r.Index)
        rec = {"pdb": r.pdb, "res1": "{}:{}{}".format(*r.key[0]), "res2": "{}:{}{}".format(*r.key[1]),
               "dssr_class": canon(r.dssr), "rubric_class": "|".join(sorted({canon(c) for c in r.rubric})),
               "votes_rubric": r.votes_rubric, "votes_dssr": r.votes_dssr}
        if residues is not None and len(r.rubric) == 1:
            bases = [names[r.pdb].get("{}|{}|{}".format(*k)) for k in r.key]
            if set(bases) <= set("ACGU"):
                rec["rmsd_dssr"] = fit(templates, r.dssr, residues, bases)
                rec["rmsd_rubric"] = fit(templates, r.rubric[0], residues, bases)
        records.append(rec)
    fit_table = pd.DataFrame(records).dropna(subset=["rmsd_dssr", "rmsd_rubric"]).reset_index(drop=True)
    fit_table["difference"] = fit_table.rmsd_dssr - fit_table.rmsd_rubric
    fit_table["closer_to"] = np.select([fit_table.difference > TIE, fit_table.difference < -TIE],
                                       ["RUBRIC", "DSSR"], "within 0.1 A")
    fit_table["conflict"] = fit_table.rubric_class + " (RUBRIC) vs " + fit_table.dssr_class + " (DSSR)"
    fit_table.to_csv(SOURCE_DATA / "template_fit_class_conflicts.csv", index=False)
    fit_table.to_csv(data_file("validation/template_fit_class_conflicts.csv"), index=False)
    print(f"{len(fit_table):,} of {len(conflict_table):,} conflicts fitted to both templates; "
          f"{(fit_table.difference == 0).sum()} exact ties")

    largest = fit_table.conflict.value_counts().index[0]
    summary = pd.DataFrame({"all conflicts": summarize(fit_table),
                            f"excluding {largest}": summarize(fit_table[fit_table.conflict != largest])}).T
    print("\nmedian RMSD to the template of RUBRIC's class:", fit_table.rmsd_rubric.median().round(2), "A")
    print("median RMSD to the template of DSSR's class:   ", fit_table.rmsd_dssr.median().round(2), "A")
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print("\n", summary.to_string(), sep="")
        by_type = (fit_table.groupby("conflict").closer_to.value_counts(normalize=True).unstack(fill_value=0)
                   .mul(100).round(1).join(fit_table.conflict.value_counts().rename("pairs")))
        print("\nby conflict type (at least 100 pairs)")
        print(by_type[by_type.pairs >= 100].sort_values("pairs", ascending=False)
              [["pairs", "RUBRIC", "DSSR", "within 0.1 A"]].to_string())
        referees = np.select([fit_table.votes_rubric > fit_table.votes_dssr,
                              fit_table.votes_rubric < fit_table.votes_dssr], ["RUBRIC", "DSSR"], "tie")
        print("\nreferee verdict (rows) against geometric verdict (columns)")
        print(pd.crosstab(pd.Series(referees, name="referees"), fit_table.closer_to, margins=True).to_string())


if __name__ == "__main__":
    main()
