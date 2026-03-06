"""Convert dot-bracket RNA secondary structure to 3D coordinates.

Uses A-form RNA helix geometry to lay out a pre-miRNA hairpin in 3D space.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import cos, pi, sin, sqrt

import numpy as np
from scipy.interpolate import CubicSpline

from viewer.config import HELIX_RADIUS, HELIX_RISE, HELIX_TWIST_DEG


@dataclass
class HairpinLayout:
    """3D layout of a pre-miRNA hairpin structure."""

    base_positions: np.ndarray  # (N, 3) xyz for each base
    base_identities: list[str]  # nucleotide letter per position
    pair_indices: list[tuple[int, int]]  # base pair (i, j) index pairs
    region_labels: list[str]  # 'stem', 'loop', 'overhang_5p', 'overhang_3p', 'bulge'
    backbone_spline: np.ndarray  # (M, 3) smooth curve through bases
    cleavage_site_positions: dict[int, np.ndarray]  # DC site -> xyz position
    sequence_length: int


def parse_dot_bracket(structure: str) -> tuple[list[tuple[int, int]], list[str]]:
    """Parse dot-bracket notation into base pairs and region labels.

    Returns:
        pairs: list of (i, j) base pair indices
        labels: per-position labels ('stem', 'loop', 'overhang_5p', 'overhang_3p')
    """
    n = len(structure)
    pairs: list[tuple[int, int]] = []
    stack: list[int] = []
    labels = [""] * n

    for i, ch in enumerate(structure):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            if stack:
                j = stack.pop()
                pairs.append((j, i))

    paired = set()
    for i, j in pairs:
        paired.add(i)
        paired.add(j)

    # Identify regions
    # Find the first and last paired positions
    if pairs:
        first_paired = min(min(i, j) for i, j in pairs)
        last_paired = max(max(i, j) for i, j in pairs)
    else:
        first_paired = 0
        last_paired = n - 1

    # Find the innermost pairs to locate the apical loop
    # Sort pairs by 5' index
    pairs_sorted = sorted(pairs, key=lambda p: p[0])
    if pairs_sorted:
        innermost = max(pairs_sorted, key=lambda p: p[0])
        loop_start = innermost[0]
        loop_end = innermost[1]
    else:
        loop_start = n // 2 - 2
        loop_end = n // 2 + 2

    for i in range(n):
        if i < first_paired:
            labels[i] = "overhang_5p"
        elif i > last_paired:
            labels[i] = "overhang_3p"
        elif i in paired:
            labels[i] = "stem"
        elif loop_start < i < loop_end:
            labels[i] = "loop"
        else:
            labels[i] = "bulge"

    return pairs, labels


def compute_hairpin_3d(
    sequence: str,
    structure: str,
    flanking_5p: int = 0,
) -> HairpinLayout:
    """Compute 3D coordinates for a pre-miRNA hairpin.

    Args:
        sequence: nucleotide sequence string
        structure: dot-bracket secondary structure
        flanking_5p: number of 5' flanking unpaired nucleotides
    """
    n = len(sequence)
    if len(structure) < n:
        structure = structure + "." * (n - len(structure))
    elif len(structure) > n:
        structure = structure[:n]

    pairs, labels = parse_dot_bracket(structure)
    positions = np.zeros((n, 3), dtype=np.float64)

    twist_rad = HELIX_TWIST_DEG * pi / 180.0

    # Build a mapping: for each base, which pair index is it?
    # Sort pairs by 5' index to get pair ordering along the stem
    pairs_sorted = sorted(pairs, key=lambda p: p[0])
    pair_map_5p: dict[int, int] = {}  # base index -> pair rank
    pair_map_3p: dict[int, int] = {}

    for rank, (i5, i3) in enumerate(pairs_sorted):
        pair_map_5p[i5] = rank
        pair_map_3p[i3] = rank

    num_pairs = len(pairs_sorted)

    # Place stem bases along helix
    for rank, (i5, i3) in enumerate(pairs_sorted):
        z = rank * HELIX_RISE
        angle_5p = rank * twist_rad
        angle_3p = angle_5p + pi  # opposite side

        positions[i5] = [
            HELIX_RADIUS * cos(angle_5p),
            HELIX_RADIUS * sin(angle_5p),
            z,
        ]
        positions[i3] = [
            HELIX_RADIUS * cos(angle_3p),
            HELIX_RADIUS * sin(angle_3p),
            z,
        ]

    # Place apical loop bases
    if pairs_sorted:
        innermost = pairs_sorted[-1]
        loop_indices = [
            i for i in range(innermost[0] + 1, innermost[1])
            if labels[i] == "loop"
        ]
        if loop_indices:
            # Arc from the 5' innermost stem base to the 3' innermost stem base
            pos_start = positions[innermost[0]]
            pos_end = positions[innermost[1]]
            top_z = (num_pairs) * HELIX_RISE + HELIX_RADIUS * 0.5

            num_loop = len(loop_indices)
            for li, idx in enumerate(loop_indices):
                t = (li + 1) / (num_loop + 1)
                # Circular arc in the plane above the helix
                mid = (pos_start + pos_end) / 2.0
                perp_x = -(pos_end[1] - pos_start[1])
                perp_y = pos_end[0] - pos_start[0]
                perp_len = sqrt(perp_x**2 + perp_y**2) or 1.0
                perp_x /= perp_len
                perp_y /= perp_len

                arc_angle = pi * t
                bulge = sin(arc_angle) * HELIX_RADIUS * 0.6

                positions[idx] = [
                    (1 - t) * pos_start[0] + t * pos_end[0] + perp_x * bulge,
                    (1 - t) * pos_start[1] + t * pos_end[1] + perp_y * bulge,
                    top_z + sin(arc_angle) * HELIX_RISE * 2,
                ]

    # Place overhang bases
    # 5' overhang: extend downward from the first stem base
    overhang_5p = [i for i in range(n) if labels[i] == "overhang_5p"]
    if overhang_5p and pairs_sorted:
        first_stem = pairs_sorted[0][0]
        base_pos = positions[first_stem].copy()
        for k, idx in enumerate(reversed(overhang_5p)):
            positions[idx] = base_pos + np.array([0, 0, -(k + 1) * HELIX_RISE])

    # 3' overhang: extend downward from the last stem 3' base
    overhang_3p = [i for i in range(n) if labels[i] == "overhang_3p"]
    if overhang_3p and pairs_sorted:
        last_stem_3p = pairs_sorted[0][1]
        base_pos = positions[last_stem_3p].copy()
        for k, idx in enumerate(overhang_3p):
            positions[idx] = base_pos + np.array([0, 0, -(k + 1) * HELIX_RISE])

    # Place bulge bases by interpolation
    bulge_indices = [i for i in range(n) if labels[i] == "bulge"]
    for idx in bulge_indices:
        # Find nearest paired neighbors
        prev_paired = idx - 1
        while prev_paired >= 0 and labels[prev_paired] not in ("stem", "loop"):
            prev_paired -= 1
        next_paired = idx + 1
        while next_paired < n and labels[next_paired] not in ("stem", "loop"):
            next_paired += 1
        if prev_paired >= 0 and next_paired < n:
            t = 0.5
            positions[idx] = (
                (1 - t) * positions[prev_paired] + t * positions[next_paired]
            )
            # Offset outward slightly for bulge visibility
            positions[idx][0] += HELIX_RADIUS * 0.3
        elif prev_paired >= 0:
            positions[idx] = positions[prev_paired] + np.array([HELIX_RISE, 0, 0])

    # Generate backbone spline (ordered 5' to 3')
    ordered_indices = list(range(n))
    if n >= 4:
        t_param = np.arange(n, dtype=np.float64)
        try:
            cs_x = CubicSpline(t_param, positions[ordered_indices, 0])
            cs_y = CubicSpline(t_param, positions[ordered_indices, 1])
            cs_z = CubicSpline(t_param, positions[ordered_indices, 2])
            t_fine = np.linspace(0, n - 1, n * 5)
            backbone = np.column_stack([cs_x(t_fine), cs_y(t_fine), cs_z(t_fine)])
        except Exception:
            backbone = positions.copy()
    else:
        backbone = positions.copy()

    # Compute cleavage site positions
    # DC_N means the product is N nt long counting from the 5' end
    # The cleavage happens between position N and N+1 (0-indexed: N-1 and N)
    cleavage_positions: dict[int, np.ndarray] = {}
    for site in [20, 21, 22, 23]:
        if site < n:
            # Position between base site-1 and base site
            idx_before = site - 1
            idx_after = min(site, n - 1)
            cleavage_positions[site] = (
                positions[idx_before] + positions[idx_after]
            ) / 2.0

    base_ids = list(sequence[:n])

    return HairpinLayout(
        base_positions=positions,
        base_identities=base_ids,
        pair_indices=pairs_sorted,
        region_labels=labels,
        backbone_spline=backbone,
        cleavage_site_positions=cleavage_positions,
        sequence_length=n,
    )
