"""Tests for viewer.chemistry.mirna_context.MiRNAContextAnalyzer."""

from __future__ import annotations

import pytest

from viewer.chemistry.mirna_context import MiRNAContextAnalyzer


class TestDominantProduct:
    """Tests for get_dominant_product."""

    def test_dominant_product_default(self):
        """Default dominant site is 21, corresponding to miR-324-5p."""
        analyzer = MiRNAContextAnalyzer()
        product = analyzer.get_dominant_product()
        assert product is not None
        assert product.name == "miR-324-5p"
        assert product.strand == "5p"
        assert product.cleavage_site == 21


class TestStrandBias:
    """Tests for get_strand_bias."""

    def test_strand_bias_5p(self):
        """dc21=0.8, dc22=0.2 reports 5p bias."""
        analyzer = MiRNAContextAnalyzer()
        bias = analyzer.get_strand_bias(dc21_acc=0.8, dc22_acc=0.2)
        assert "5p" in bias.lower()

    def test_strand_bias_balanced(self):
        """dc21 ~= dc22 reports balanced."""
        analyzer = MiRNAContextAnalyzer()
        bias = analyzer.get_strand_bias(dc21_acc=0.50, dc22_acc=0.50)
        assert "balanced" in bias.lower()


class TestAffectedGenes:
    """Tests for get_affected_genes."""

    def test_affected_genes(self):
        """Default dominant product (miR-324-5p) has target genes."""
        analyzer = MiRNAContextAnalyzer()
        genes = analyzer.get_affected_genes()
        assert len(genes) > 0
        gene_symbols = [g.symbol for g in genes]
        assert "GLI1" in gene_symbols


class TestGeneImpact:
    """Tests for predict_modification_gene_impact."""

    def test_gene_impact_minimal(self):
        """dc_ratio_shift=0.0 produces 'minimal' impact descriptions."""
        analyzer = MiRNAContextAnalyzer()
        impacts = analyzer.predict_modification_gene_impact(dc_ratio_shift=0.0)
        assert len(impacts) > 0
        for gene, description in impacts:
            assert "minimal" in description.lower()
