"""Policy for which H-bond distances define a base pair's distance score."""

from __future__ import annotations

from collections.abc import Iterable

MAX_STRONG_DISTANCE = 3.5

MIN_PHYSICAL_DISTANCE = 2.2

BASE_ATOMS: frozenset[str] = frozenset({
    "N1", "N2", "N3", "N4", "N6", "N7", "N9", "O2", "O4", "O6",
})


def distances_to_score(distances: Iterable[float | None],
                       bp_type: str | None = None, lw_class: str | None = None) -> list[float]:
    """Every base-base bond of at least MIN_PHYSICAL_DISTANCE, or none when the pair has
    no bond within MAX_STRONG_DISTANCE (the same rule for every class)."""
    phys = [d for d in distances if d is not None and d >= MIN_PHYSICAL_DISTANCE]
    if not any(d <= MAX_STRONG_DISTANCE for d in phys):
        return []
    return phys


def select_distance_rows(hb_df):
    """The H-bond rows the distance tables are built from, chosen as distances_to_score does."""
    df = hb_df[
        hb_df["donor_atom"].isin(BASE_ATOMS)
        & hb_df["acceptor_atom"].isin(BASE_ATOMS)
        & (hb_df["distance"] >= MIN_PHYSICAL_DISTANCE)
    ]
    has_strong = (df["distance"] <= MAX_STRONG_DISTANCE).groupby(
        [df["pdb_id"], df["pair_idx"]]).transform("any")
    return df[has_strong]


def base_base_rows(hb_df):
    return hb_df[
        hb_df["donor_atom"].isin(BASE_ATOMS)
        & hb_df["acceptor_atom"].isin(BASE_ATOMS)
    ]
