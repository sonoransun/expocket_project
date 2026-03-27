"""Configuration constants, color schemes, and default paths."""

from __future__ import annotations

# Nucleotide colors (from existing analysis scripts)
NT_COLORS = {
    "A": "#89CFF0",
    "T": "#3da4dc",
    "U": "#3da4dc",
    "G": "#9B111E",
    "C": "#5AC9A1",
}

NT_COLORS_FLOAT = {
    "A": (0.537, 0.812, 0.941, 1.0),
    "T": (0.239, 0.643, 0.863, 1.0),
    "U": (0.239, 0.643, 0.863, 1.0),
    "G": (0.608, 0.067, 0.118, 1.0),
    "C": (0.353, 0.788, 0.631, 1.0),
}

CLEAVAGE_COLORS = {
    20: "#ca5cdd",
    21: "#5AC9A1",
    22: "#f94449",
    23: "#FFAE42",
}

CLEAVAGE_COLORS_FLOAT = {
    20: (0.792, 0.361, 0.867, 0.6),
    21: (0.353, 0.788, 0.631, 0.6),
    22: (0.976, 0.267, 0.286, 0.6),
    23: (1.0, 0.682, 0.259, 0.6),
}

ENZYME_COLORS = {
    "hdicer": "#5AC9A1",
    "dcr": "#FFAE42",
}

# Groups (5' terminal nucleotide)
GROUPS = ["A", "T", "G", "C"]
CLEAVAGE_SITES = [20, 21, 22, 23]

# A-form RNA helix geometry parameters (angstroms / degrees)
HELIX_RISE = 2.81
HELIX_TWIST_DEG = 32.7
HELIX_RADIUS = 9.4

# Default data paths (from existing scripts)
DEFAULT_DATA_PATHS = {
    "human_folder": "D:/1.miRNA/END_randomized project/dataframes",
    "fly_folder": "D:/1.miRNA/END_randomized project/DCR1-DCL1/rawcount_dcr_pnk/",
    "structure_file": "Pre-mir-324-end-randomization-structure.bed",
}

# Canonical pre-mir-324 reference (for mock data)
PREMIR324_REFERENCE_SEQ = (
    "CGCAUCCCCUAGGGCAUUGGUGUGAAAGCUGGAGCCUCUGAAUCCUUGCUUUC"
    "CCCUAGCGUCG"
)
PREMIR324_REFERENCE_STRUCT = (
    ".((((((((((((((((((((....))))))))))).)))))))))."
)

# Property colormaps for 2D/3D visualization
PROPERTY_COLORMAPS = {
    "molecular_weight": "YlOrRd",
    "vdw_volume": "YlOrRd",
    "h_bond_donors": "Blues",
    "h_bond_acceptors": "Greens",
    "pka": "coolwarm",
    "charge_at_ph7": "RdBu_r",
    "hydrophobicity_index": "RdYlGn_r",
    "sugar_pucker_preference": "PuBuGn",
    "stacking_energy_5p": "viridis",
    "stacking_energy_3p": "viridis",
    "propeller_twist": "coolwarm",
    "buckle": "coolwarm",
}

# DICER pocket domain colors (RGBA float)
DOMAIN_COLORS = {
    "PAZ": (0.2, 0.6, 0.9, 0.8),
    "platform": (0.5, 0.8, 0.3, 0.8),
    "RNaseIIIa": (0.9, 0.2, 0.2, 0.8),
    "RNaseIIIb": (0.9, 0.5, 0.1, 0.8),
    "helicase": (0.7, 0.4, 0.9, 0.8),
}

# Modification display colors (RGBA float, from modification_db)
MODIFICATION_COLORS = {
    "2OMe": (0.2, 0.8, 0.9, 0.85),
    "LNA": (0.9, 0.3, 0.1, 0.85),
    "PSI": (0.7, 0.5, 0.9, 0.85),
    "m6A": (0.9, 0.8, 0.2, 0.85),
    "2F": (0.3, 0.9, 0.4, 0.85),
    "PS": (1.0, 0.6, 0.0, 0.85),
    "s4U": (0.8, 0.7, 0.1, 0.85),
    "m5C": (0.5, 0.8, 0.7, 0.85),
    "INO": (0.6, 0.4, 0.8, 0.85),
}

# Synthesis yield display thresholds and colors
YIELD_COLORS = {
    "good": "#5AC9A1",   # > 50% yield
    "warn": "#FFAE42",   # 20-50% yield
    "poor": "#f44",      # < 20% yield
}
YIELD_THRESHOLDS = {"good": 0.5, "warn": 0.2}

# Gene editing tool colors (RGBA float)
EDITING_TOOL_COLORS: dict[str, tuple[float, float, float, float]] = {
    "cas9":      (0.2,  0.6,  1.0,  0.85),  # blue
    "cas12a":    (0.4,  0.9,  0.4,  0.85),  # green
    "casclover": (1.0,  0.7,  0.1,  0.85),  # amber
    "talen":     (0.9,  0.3,  0.3,  0.85),  # red
    "zfn":       (0.8,  0.4,  0.9,  0.85),  # purple
    "retron":    (0.3,  0.9,  0.8,  0.85),  # teal
    "nicer":     (1.0,  0.5,  0.2,  0.85),  # orange
}

# Synergy display colors for double replacement screening
SYNERGY_COLORS = {
    "cooperative": "#5AC9A1",    # synergy > 0
    "antagonistic": "#f44",      # synergy < 0
    "neutral": "#888",           # synergy ~ 0
}
