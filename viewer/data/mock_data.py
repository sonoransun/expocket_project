"""Generate realistic mock data for development and demo purposes."""

from __future__ import annotations

import itertools
import math
import random

import numpy as np

from viewer.config import GROUPS, PREMIR324_REFERENCE_SEQ, PREMIR324_REFERENCE_STRUCT
from viewer.data.schema import CleavageRecord, VariantDataset, VariantInfo

_NTS = ["A", "T", "G", "C"]

# Biological priors: group-dependent DC21 accuracy tendencies
# A-group tends higher DC21, G-group tends higher DC22 (from hsa_324_data_khoa.py)
_GROUP_DC21_MEAN = {"A": 0.55, "T": 0.50, "G": 0.40, "C": 0.45}
_GROUP_DC22_MEAN = {"A": 0.20, "T": 0.25, "G": 0.35, "C": 0.30}


def _make_sequence(group: str, rand_nts: str) -> str:
    """Build a pre-miRNA sequence from the reference, substituting terminal nts."""
    seq = list(PREMIR324_REFERENCE_SEQ)
    # Replace 5' nucleotide
    seq[0] = group
    # Replace last 3 nucleotides (3' randomized end)
    for i, nt in enumerate(rand_nts):
        seq[-(3 - i)] = nt
    return "".join(seq)


def _make_structure(group: str, rand_nts: str) -> str:
    """Generate a dot-bracket structure with slight variation."""
    base = PREMIR324_REFERENCE_STRUCT
    # Most variants share the same structure; ~30% get a 1-nt 5' overhang variant
    if random.random() < 0.3:
        # Convert leading paired base to unpaired
        base = "." + base[1:]
    return base


def _generate_accuracy(group: str, rand_nts: str, site: int) -> float:
    """Generate a plausible cleavage accuracy value."""
    if site == 21:
        mean = _GROUP_DC21_MEAN[group]
        # 3' opposite nucleotide (first of randomized) modulates accuracy
        opp = rand_nts[0]
        if group == "A" and opp == "T":
            mean += 0.10
        elif group == "G" and opp == "C":
            mean += 0.08
    elif site == 22:
        mean = _GROUP_DC22_MEAN[group]
        opp = rand_nts[0]
        if group == "G" and opp == "C":
            mean += 0.10
    elif site == 20:
        mean = 0.08
    else:  # 23
        mean = 0.06

    val = random.gauss(mean, 0.12)
    return max(0.0, min(1.0, val))


def generate_mock_dataset(enzyme: str = "hdicer", seed: int = 42) -> VariantDataset:
    """Generate a complete mock dataset with 256 variants."""
    random.seed(seed)
    np.random.seed(seed)

    variants: list[VariantInfo] = []
    cleavage_data: dict[str, list[CleavageRecord]] = {}

    for group in GROUPS:
        for combo in itertools.product(_NTS, repeat=3):
            rand_nts = "".join(combo)
            variant_id = f"{group}_{rand_nts}"

            seq = _make_sequence(group, rand_nts)
            struct = _make_structure(group, rand_nts)
            flank = 1 if struct[0] == "." else 0

            vi = VariantInfo(
                variant=variant_id,
                group=group,
                randomized_nts=rand_nts,
                pre_mirna_sequence=seq,
                concrete_struct=struct,
                flanking_length_5p=flank,
            )
            variants.append(vi)

            records = []
            accuracies = {}
            for site in [20, 21, 22, 23]:
                acc = _generate_accuracy(group, rand_nts, site)
                accuracies[site] = acc

            # Normalize so accuracies sum to ~1.0
            total = sum(accuracies.values())
            if total > 0:
                for site in accuracies:
                    accuracies[site] /= total

            for site in [20, 21, 22, 23]:
                acc = accuracies[site]
                # Simulate 3 replicates with noise
                reps = [max(0.0, min(1.0, acc + random.gauss(0, 0.03))) for _ in range(3)]
                eff = math.log2(acc + 0.1) - math.log2(0.5)
                records.append(
                    CleavageRecord(
                        variant=variant_id,
                        cleavage_site=site,
                        mean_accuracy=acc,
                        mean_positional_efficiency=eff,
                        mean_global_efficiency=eff * 0.8,
                        accuracy_rep1=reps[0],
                        accuracy_rep2=reps[1],
                        accuracy_rep3=reps[2],
                    )
                )
            cleavage_data[variant_id] = records

    return VariantDataset(variants=variants, cleavage_data=cleavage_data, enzyme=enzyme)
