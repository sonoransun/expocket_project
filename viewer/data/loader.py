"""Load experimental data from TSV files into VariantDataset."""

from __future__ import annotations

import logging
import os
from pathlib import Path

import pandas as pd

from viewer.config import DEFAULT_DATA_PATHS, GROUPS
from viewer.data.mock_data import generate_mock_dataset
from viewer.data.schema import CleavageRecord, VariantDataset, VariantInfo

_log = logging.getLogger(__name__)


def _dataframe_to_dataset(df: pd.DataFrame, enzyme: str) -> VariantDataset:
    """Convert a combined pandas DataFrame to a VariantDataset."""
    if "Variant" not in df.columns:
        _log.warning("Missing required 'Variant' column — returning empty dataset")
        return VariantDataset(variants=[], cleavage_data={}, enzyme=enzyme)

    variants: list[VariantInfo] = []
    cleavage_data: dict[str, list[CleavageRecord]] = {}

    # Extract unique variants
    variant_cols = [
        "Variant", "Group", "Randomized_nts", "Pre_miRNA_sequence",
        "concrete_struct", "5p_flanking_length",
    ]
    opt_cols = ["New_define_structure_1", "New_define_structure_2"]
    available = [c for c in variant_cols if c in df.columns]

    seen_variants: set[str] = set()
    skipped = 0
    for row_idx, row in df.iterrows():
        try:
            vid = row["Variant"]
            if vid not in seen_variants:
                seen_variants.add(vid)
                vi = VariantInfo(
                    variant=vid,
                    group=row.get("Group", ""),
                    randomized_nts=row.get("Randomized_nts", ""),
                    pre_mirna_sequence=row.get("Pre_miRNA_sequence", ""),
                    concrete_struct=row.get("concrete_struct", ""),
                    flanking_length_5p=int(row.get("5p_flanking_length", 0)),
                    new_define_structure_1=row.get("New_define_structure_1", ""),
                    new_define_structure_2=row.get("New_define_structure_2", ""),
                )
                variants.append(vi)

            # Cleavage record
            if "Cleavage_site" in df.columns:
                rec = CleavageRecord(
                    variant=vid,
                    cleavage_site=int(row.get("Cleavage_site", 0)),
                    mean_accuracy=float(row.get("Mean_Cleavage_accuracy", 0)),
                    mean_positional_efficiency=float(
                        row.get("Mean_Position_efficiency", 0)
                    ),
                    mean_global_efficiency=float(
                        row.get("Mean_Global_efficiency", 0)
                    ),
                    accuracy_rep1=float(row.get("Cleavage_accuracy_rep1", 0)),
                    accuracy_rep2=float(row.get("Cleavage_accuracy_rep2", 0)),
                    accuracy_rep3=float(row.get("Cleavage_accuracy_rep3", 0)),
                )
                cleavage_data.setdefault(vid, []).append(rec)
        except (KeyError, ValueError, TypeError) as exc:
            skipped += 1
            _log.warning("Skipping malformed row %s: %s", row_idx, exc)
            continue

    if skipped:
        _log.info("Skipped %d malformed rows out of %d total", skipped, len(df))

    return VariantDataset(variants=variants, cleavage_data=cleavage_data, enzyme=enzyme)


def load_human_data(folder_path: str | None = None) -> VariantDataset:
    """Load human DICER PNK data from 4 group TSV files.

    Falls back to mock data if files are not found.
    """
    folder = folder_path or DEFAULT_DATA_PATHS["human_folder"]
    dfs = []
    for group in GROUPS:
        path = os.path.join(folder, f"human_df{group}_pnk.txt")
        if not os.path.exists(path):
            return generate_mock_dataset(enzyme="hdicer")
        try:
            dfs.append(pd.read_csv(path, sep="\t"))
        except (pd.errors.ParserError, pd.errors.EmptyDataError, UnicodeDecodeError) as exc:
            _log.warning("Failed to parse %s: %s — falling back to mock data", path, exc)
            return generate_mock_dataset(enzyme="hdicer")

    combined = pd.concat(dfs, ignore_index=True)
    return _dataframe_to_dataset(combined, enzyme="hdicer")


def load_fly_data(folder_path: str | None = None) -> VariantDataset:
    """Load fly DCR-1 PNK combined data.

    Falls back to mock data if the file is not found.
    """
    folder = folder_path or DEFAULT_DATA_PATHS["fly_folder"]
    path = os.path.join(folder, "df_dcr_pnk_combine.txt")
    if not os.path.exists(path):
        return generate_mock_dataset(enzyme="dcr", seed=99)

    try:
        df = pd.read_csv(path, sep="\t")
    except (pd.errors.ParserError, pd.errors.EmptyDataError, UnicodeDecodeError) as exc:
        _log.warning("Failed to parse %s: %s — falling back to mock data", path, exc)
        return generate_mock_dataset(enzyme="dcr", seed=99)
    return _dataframe_to_dataset(df, enzyme="dcr")


def load_dataset(
    folder_path: str | None = None, enzyme: str = "hdicer"
) -> VariantDataset:
    """Load a dataset by enzyme type, falling back to mock data."""
    if enzyme == "hdicer":
        return load_human_data(folder_path)
    elif enzyme == "dcr":
        return load_fly_data(folder_path)
    else:
        return generate_mock_dataset(enzyme=enzyme)
