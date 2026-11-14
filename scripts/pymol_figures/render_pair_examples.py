"""Render the base-pair example panels of Fig. S8 and Fig. S9 with PyMOL.

Each pair is drawn face-on to its mean base plane, nucleotide 1 to 2 along the screen x
axis, at the camera distance of the Figure 4 pair sessions, so every panel prints at the
same base size as Figures 3 and 4. For every example the script writes

  data/pymol/renders/<key>.png        the render (2100 x 1575 px)
  data/pymol/renders/<key>_meta.json  H-bond lengths and pixel anchors for their labels
  data/pymol/sessions/<key>.pse       the PyMOL session

The examples are read from the EXAMPLES table of each figure's notebook, so the notebook
stays the single place where the pairs are defined. The notebooks only read the files
above, so they run without PyMOL or the structure files. Fig. S8 shows bases and riboses
(thin ribose bonds); Fig. S9 shows bases with the C1' stub only, as in Figures 3 and 4.

Needs PyMOL and the mmCIF files in data/structures (see scripts/fetch_structures.py).

Usage (from the repo root):  python scripts/pymol_figures/render_pair_examples.py [S8] [S9]
"""
import ast
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image

try:
    from pymol import cmd
except ImportError:  # importable without PyMOL (CI imports every script); main() needs it
    cmd = None

ROOT = Path(__file__).resolve().parents[2]
STRUCTURES = ROOT / "data" / "structures"
RENDERS = ROOT / "data" / "pymol" / "renders"
SESSIONS = ROOT / "data" / "pymol" / "sessions"
NOTEBOOKS = ROOT / "notebooks" / "supplemental_figures"

W, H = 2100, 1575
CAMERA_DIST = 75.24  # camera distance of the Figure 4 pair sessions (6D_orig.pse)
RIBOSE = "C1'+C2'+C3'+C4'+O4'+O2'"
BASE_COLORS = {"G": [0.878, 0.0, 0.0], "C": [0.0, 0.565, 0.125],     # Figures 3B and 4K/L
               "A": [0.941, 0.502, 0.0], "U": [0.0, 0.0, 0.878]}
RIBOSE_RADIUS = 0.10

# figure: (notebook, atoms hidden, thin ribose bonds)
FIGURES = {
    "S8": ("figS8.ipynb", "P+OP1+OP2+OP3+O5'+C5'+O3'", True),
    "S9": ("figS9.ipynb", "P+OP1+OP2+OP3+O5'+C5'+O3'+C2'+C3'+C4'+O4'+O2'", False),
}


def examples(notebook):
    """The literal EXAMPLES dict of a figure notebook."""
    nb = json.loads((NOTEBOOKS / notebook).read_text())
    for cell in nb["cells"]:
        if cell["cell_type"] != "code":
            continue
        for node in ast.parse("".join(cell["source"])).body:
            if isinstance(node, ast.Assign) and any(getattr(t, "id", None) == "EXAMPLES" for t in node.targets):
                return ast.literal_eval(node.value)
    raise SystemExit(f"no EXAMPLES table in {notebook}")


def _sel(nt, hide):
    ch, resi, resn = nt
    return f"s and chain {ch} and resi {resi} and resn {resn} and not name {hide} and not elem H"


def _base_normal(obj):
    xyz = np.array(cmd.get_coords(f"{obj} and not name {RIBOSE}"))
    return np.linalg.svd(xyz - xyz.mean(axis=0))[2][2]


def _setup():
    cmd.set("bg_rgb", "white"); cmd.set("ray_opaque_background", 1)
    cmd.set("stick_radius", 0.18); cmd.set("valence", 1); cmd.set("valence_size", 0.06)
    cmd.set("antialias", 2); cmd.set("ray_shadows", 0); cmd.set("specular", 1)
    cmd.set("ambient", 0.14); cmd.set("direct", 0.45); cmd.set("reflect", 0.45)
    cmd.set("depth_cue", 1); cmd.set("field_of_view", 20); cmd.set("orthoscopic", 0)
    cmd.set("dash_gap", 0.35); cmd.set("dash_length", 0.15); cmd.set("dash_width", 2.5)
    cmd.set("dash_round_ends", 1); cmd.set("dash_radius", 0.0)
    for name, rgb in BASE_COLORS.items():
        cmd.set_color(f"base_{name}", rgb)


def _colour():
    for resn in BASE_COLORS:
        cmd.color(f"base_{resn}", f"(n1 or n2) and resn {resn}")
    cmd.color("black", "d*")


def _locate(obj, key):
    """Pixel centroid and lowest pixel of one object, rendered alone in pure blue."""
    tmp = RENDERS / f"_{key}_locate.png"
    shown = cmd.get_names("public_objects", enabled_only=1)
    cmd.disable("all"); cmd.enable(obj)
    cmd.color("blue", obj)
    cmd.png(str(tmp), width=W, height=H, ray=1)
    _colour()
    for o in shown:
        cmd.enable(o)
    a = np.asarray(Image.open(tmp).convert("RGB")).astype(int)
    tmp.unlink()
    ys, xs = np.nonzero((a[..., 2] > 150) & (a[..., 0] < 90) & (a[..., 1] < 90))
    return [float(xs.mean()), float(ys.mean()), float(ys.max())]


def render(key, example, hide, thin_ribose):
    pdb, nt1, nt2, dashes, _, _ = example
    cmd.reinitialize(); _setup()
    cmd.load(str(STRUCTURES / f"{pdb}.cif"), "s")
    cmd.create("n1", _sel(nt1, hide)); cmd.create("n2", _sel(nt2, hide)); cmd.delete("s")
    cmd.hide("everything"); cmd.show("sticks", "n1 or n2")
    if thin_ribose:
        cmd.set_bond("stick_radius", RIBOSE_RADIUS, f"(n1 or n2) and name {RIBOSE}")
    _colour()

    # face-on to the mean base plane, nucleotide 1 -> 2 along screen x
    c1, c2 = (np.array(cmd.get_coords(o)).mean(axis=0) for o in ("n1", "n2"))
    n1, n2 = _base_normal("n1"), _base_normal("n2")
    z = n1 + (n2 if n1 @ n2 >= 0 else -n2); z /= np.linalg.norm(z)
    x = (c2 - c1) - ((c2 - c1) @ z) * z; x /= np.linalg.norm(x)
    y = np.cross(z, x)
    centre = (c1 + c2) / 2
    rot = np.column_stack([x, y, z])
    cmd.set_view(list(rot.flatten()) + [0.0, 0.0, -CAMERA_DIST] + list(centre)
                 + [CAMERA_DIST - 20.0, CAMERA_DIST + 20.0, -20.0])

    meta = {"dist": {}, "px": {}}
    for i, (a, b) in enumerate(dashes):
        name = f"d{i}"
        meta["dist"][name] = round(cmd.distance(name, f"n1 and name {a}", f"n2 and name {b}"), 2)
        cmd.hide("labels", name)
    _colour()
    for obj in ["n1", "n2"] + [f"d{i}" for i in range(len(dashes))]:
        meta["px"][obj] = _locate(obj, key)
    cmd.png(str(RENDERS / f"{key}.png"), width=W, height=H, dpi=300, ray=1)
    (RENDERS / f"{key}_meta.json").write_text(json.dumps(meta, indent=1) + "\n")
    cmd.save(str(SESSIONS / f"{key}.pse"))
    return meta["dist"]


def main():
    if cmd is None:
        raise SystemExit("PyMOL is required: install it (e.g. conda install -c conda-forge pymol-open-source)")
    wanted = [a.upper() for a in sys.argv[1:]] or list(FIGURES)
    unknown = [w for w in wanted if w not in FIGURES]
    if unknown:
        raise SystemExit(f"unknown figure(s) {unknown}; choose from {list(FIGURES)}")
    RENDERS.mkdir(parents=True, exist_ok=True)
    SESSIONS.mkdir(parents=True, exist_ok=True)
    for fig in wanted:
        notebook, hide, thin_ribose = FIGURES[fig]
        for key, example in examples(notebook).items():
            print(fig, key, render(key, example, hide, thin_ribose))


if __name__ == "__main__":
    main()
