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
