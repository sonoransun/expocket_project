"""Shared test fixtures for expocket viewer tests."""

from __future__ import annotations

import pytest
import numpy as np

from viewer.config import PREMIR324_REFERENCE_SEQ, PREMIR324_REFERENCE_STRUCT
from viewer.data.schema import (
    CleavageRecord,
    EnrichedVariantDataset,
    ModificationState,
    VariantDataset,
    VariantInfo,
)


def _make_variant(group: str, nts: str) -> VariantInfo:
    """Build a variant with known sequence from reference."""
    seq = list(PREMIR324_REFERENCE_SEQ)
    seq[0] = group
    for i, nt in enumerate(nts):
        seq[-(3 - i)] = nt
    return VariantInfo(
        variant=f"{group}_{nts}",
        group=group,
        randomized_nts=nts,
        pre_mirna_sequence="".join(seq),
        concrete_struct=PREMIR324_REFERENCE_STRUCT,
        flanking_length_5p=0,
    )


def _make_cleavage(variant_id: str, acc20: float, acc21: float,
                   acc22: float, acc23: float) -> list[CleavageRecord]:
    return [
        CleavageRecord(variant=variant_id, cleavage_site=20, mean_accuracy=acc20),
        CleavageRecord(variant=variant_id, cleavage_site=21, mean_accuracy=acc21),
        CleavageRecord(variant=variant_id, cleavage_site=22, mean_accuracy=acc22),
        CleavageRecord(variant=variant_id, cleavage_site=23, mean_accuracy=acc23),
    ]


@pytest.fixture
def sample_variant_info() -> VariantInfo:
    return _make_variant("A", "AAA")


@pytest.fixture
def sample_cleavage_records() -> list[CleavageRecord]:
    return _make_cleavage("A_AAA", 0.10, 0.55, 0.25, 0.10)


@pytest.fixture
def small_dataset() -> VariantDataset:
    """4 variants, one per group, with cleavage data."""
    variants = [
        _make_variant("A", "AAA"),
        _make_variant("T", "TTT"),
        _make_variant("G", "GGG"),
        _make_variant("C", "CCC"),
    ]
    cleavage_data = {}
    cleavage_data["A_AAA"] = _make_cleavage("A_AAA", 0.10, 0.55, 0.25, 0.10)
    cleavage_data["T_TTT"] = _make_cleavage("T_TTT", 0.08, 0.50, 0.30, 0.12)
    cleavage_data["G_GGG"] = _make_cleavage("G_GGG", 0.12, 0.40, 0.35, 0.13)
    cleavage_data["C_CCC"] = _make_cleavage("C_CCC", 0.09, 0.45, 0.32, 0.14)
    return VariantDataset(variants=variants, cleavage_data=cleavage_data, enzyme="hdicer")


@pytest.fixture(scope="session")
def enriched_dataset() -> EnrichedVariantDataset:
    """4-variant enriched dataset (session-scoped for speed)."""
    from viewer.data.mock_data import enrich_dataset
    variants = [
        _make_variant("A", "AAA"),
        _make_variant("T", "TTT"),
        _make_variant("G", "GGG"),
        _make_variant("C", "CCC"),
    ]
    cleavage_data = {}
    cleavage_data["A_AAA"] = _make_cleavage("A_AAA", 0.10, 0.55, 0.25, 0.10)
    cleavage_data["T_TTT"] = _make_cleavage("T_TTT", 0.08, 0.50, 0.30, 0.12)
    cleavage_data["G_GGG"] = _make_cleavage("G_GGG", 0.12, 0.40, 0.35, 0.13)
    cleavage_data["C_CCC"] = _make_cleavage("C_CCC", 0.09, 0.45, 0.32, 0.14)
    ds = VariantDataset(variants=variants, cleavage_data=cleavage_data, enzyme="hdicer")
    return enrich_dataset(ds)


@pytest.fixture(scope="session")
def mock_dataset_256() -> EnrichedVariantDataset:
    """Full 256-variant enriched dataset (session-scoped)."""
    from viewer.data.mock_data import generate_mock_dataset, enrich_dataset
    ds = generate_mock_dataset(enzyme="hdicer", seed=42)
    return enrich_dataset(ds)


@pytest.fixture
def empty_dataset() -> VariantDataset:
    return VariantDataset(variants=[], cleavage_data={}, enzyme="hdicer")
