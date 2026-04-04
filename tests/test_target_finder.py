"""Tests for viewer.genome_editing.target_finder."""

from __future__ import annotations

from viewer.genome_editing.target_finder import TargetFinder
from viewer.genome_editing.tools import TargetSite


def test_find_sites_cas9(enriched_dataset):
    finder = TargetFinder()
    sites = finder.find_sites("A_AAA", "cas9", enriched_dataset)
    assert isinstance(sites, list)
    assert len(sites) > 0
    assert all(isinstance(s, TargetSite) for s in sites)
    assert all(s.tool == "cas9" for s in sites)


def test_find_sites_unknown_variant(enriched_dataset):
    finder = TargetFinder()
    sites = finder.find_sites("NONEXISTENT_999", "cas9", enriched_dataset)
    assert sites == []


def test_find_all_tools_keys(enriched_dataset):
    finder = TargetFinder()
    result = finder.find_all_tools("A_AAA", enriched_dataset)
    assert isinstance(result, dict)
    expected_keys = {"cas9", "cas12a", "casclover", "talen", "zfn", "retron", "nicer"}
    assert set(result.keys()) == expected_keys


def test_sites_sorted_by_score(enriched_dataset):
    finder = TargetFinder()
    sites = finder.find_sites("A_AAA", "cas9", enriched_dataset)
    if len(sites) >= 2:
        for i in range(len(sites) - 1):
            assert sites[i].binding_score >= sites[i + 1].binding_score


def test_target_site_has_required_fields(enriched_dataset):
    finder = TargetFinder()
    sites = finder.find_sites("A_AAA", "cas9", enriched_dataset)
    assert len(sites) > 0
    site = sites[0]
    assert hasattr(site, "tool")
    assert hasattr(site, "strand")
    assert hasattr(site, "spacer_start")
    assert hasattr(site, "binding_score")
    assert site.tool == "cas9"
    assert site.strand in ("sense", "antisense")
    assert isinstance(site.spacer_start, int)
    assert 0.0 <= site.binding_score <= 1.0
