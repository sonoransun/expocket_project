"""miRNA product context: strand bias, target genes, expression impact."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class MiRNAProduct:
    """A mature miRNA product from DICER cleavage."""

    name: str
    sequence: str
    strand: str  # "5p" or "3p"
    cleavage_site: int  # DC20, DC21, DC22, DC23


@dataclass(frozen=True)
class TargetGene:
    """A predicted miRNA target gene."""

    symbol: str
    full_name: str
    binding_type: str  # "seed", "3'-supplementary", "centered"
    context_score: float  # TargetScan-like score (-1 to 0, more negative = stronger)
    pathway: str


# Mock miR-324 products for DC21 and DC22
MIRNA_PRODUCTS = {
    21: MiRNAProduct(
        name="miR-324-5p",
        sequence="CGCAUCCCCUAGGGCAUUGGUG",
        strand="5p",
        cleavage_site=21,
    ),
    22: MiRNAProduct(
        name="miR-324-3p",
        sequence="CCAGUAUCUGAAUCCUUGCUUUC",
        strand="3p",
        cleavage_site=22,
    ),
}

# Mock target genes for miR-324 products
TARGET_GENES: dict[str, list[TargetGene]] = {
    "miR-324-5p": [
        TargetGene("GLI1", "GLI family zinc finger 1", "seed", -0.45, "Hedgehog"),
        TargetGene("SMO", "Smoothened", "seed", -0.38, "Hedgehog"),
        TargetGene("NOTCH1", "Notch receptor 1", "3'-supplementary", -0.28, "Notch"),
        TargetGene("WNT2B", "Wnt family member 2B", "seed", -0.22, "Wnt"),
        TargetGene("CCND1", "Cyclin D1", "centered", -0.18, "Cell cycle"),
    ],
    "miR-324-3p": [
        TargetGene("PTCH1", "Patched 1", "seed", -0.52, "Hedgehog"),
        TargetGene("SUFU", "SUFU negative regulator of Hedgehog", "seed", -0.35, "Hedgehog"),
        TargetGene("DVL2", "Dishevelled segment polarity protein 2", "3'-supplementary", -0.25, "Wnt"),
        TargetGene("RBPJ", "Recombination signal binding protein for Ig kappa J", "seed", -0.20, "Notch"),
    ],
}


@dataclass
class MiRNAContextAnalyzer:
    """Analyzes miRNA products and their biological context."""

    _dominant_site: int = 21

    def set_dominant_site(self, site: int) -> None:
        self._dominant_site = site

    def get_dominant_product(self) -> MiRNAProduct | None:
        return MIRNA_PRODUCTS.get(self._dominant_site)

    def get_strand_bias(self, dc21_acc: float, dc22_acc: float) -> str:
        """Return strand bias description."""
        ratio = dc21_acc / (dc22_acc + 1e-6)
        if ratio > 1.5:
            return "Strong 5p bias (miR-324-5p dominant)"
        elif ratio > 1.1:
            return "Moderate 5p bias"
        elif ratio < 0.67:
            return "Strong 3p bias (miR-324-3p dominant)"
        elif ratio < 0.9:
            return "Moderate 3p bias"
        return "Balanced 5p/3p"

    def get_affected_genes(self) -> list[TargetGene]:
        """Return target genes for the dominant product."""
        product = self.get_dominant_product()
        if product is None:
            return []
        return TARGET_GENES.get(product.name, [])

    def predict_modification_gene_impact(
        self, dc_ratio_shift: float
    ) -> list[tuple[TargetGene, str]]:
        """Predict how a modification shifts gene targeting.

        Returns (gene, impact_description) pairs.
        """
        results: list[tuple[TargetGene, str]] = []

        # If ratio shifts toward 5p, 5p targets upregulated (more suppression)
        genes_5p = TARGET_GENES.get("miR-324-5p", [])
        genes_3p = TARGET_GENES.get("miR-324-3p", [])

        if abs(dc_ratio_shift) < 0.01:
            return [(g, "minimal change") for g in genes_5p + genes_3p]

        for gene in genes_5p:
            if dc_ratio_shift > 0:
                impact = f"increased suppression (shift +{dc_ratio_shift:.3f})"
            else:
                impact = f"decreased suppression (shift {dc_ratio_shift:.3f})"
            results.append((gene, impact))

        for gene in genes_3p:
            if dc_ratio_shift < 0:
                impact = f"increased suppression (shift {dc_ratio_shift:.3f})"
            else:
                impact = f"decreased suppression (shift +{dc_ratio_shift:.3f})"
            results.append((gene, impact))

        return results
