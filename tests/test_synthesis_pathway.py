"""Tests for viewer.chemistry.synthesis_pathway.SynthesisPlanner."""

from __future__ import annotations

import pytest

from viewer.chemistry.synthesis_pathway import SynthesisPlanner


class TestPlanUnmodified:
    """Tests for planning synthesis of unmodified variants."""

    def test_plan_unmodified(self, enriched_dataset):
        """All steps use standard monomers; total_yield close to 0.995^seq_len."""
        planner = SynthesisPlanner(enriched_dataset)
        plan = planner.plan_synthesis("A_AAA", modifications={})
        seq_len = len(plan.sequence)

        # Every step should use an unmodified monomer
        for step in plan.steps:
            assert step.modification is None

        # Standard coupling efficiency is 0.995
        expected_yield = 0.995 ** seq_len
        assert abs(plan.total_yield - expected_yield) < 1e-6

    def test_plan_unknown_variant(self, enriched_dataset):
        """Unknown variant returns a plan with empty sequence and no steps."""
        planner = SynthesisPlanner(enriched_dataset)
        plan = planner.plan_synthesis("NONEXISTENT")
        assert plan.sequence == ""
        assert plan.steps == []


class TestPlanStepDetails:
    """Tests for step count, ordering, and cumulative yield."""

    def test_plan_step_count(self, enriched_dataset):
        """Number of steps equals the sequence length."""
        planner = SynthesisPlanner(enriched_dataset)
        plan = planner.plan_synthesis("A_AAA", modifications={})
        seq_len = len(plan.sequence)
        assert len(plan.steps) == seq_len

    def test_plan_3p_to_5p_order(self, enriched_dataset):
        """First step has last position (3' end), last step has position 0 (5' end)."""
        planner = SynthesisPlanner(enriched_dataset)
        plan = planner.plan_synthesis("A_AAA", modifications={})
        seq_len = len(plan.sequence)
        # Synthesis goes 3'->5', so first step is last sequence position
        assert plan.steps[0].position == seq_len - 1
        assert plan.steps[-1].position == 0

    def test_plan_cumulative_yield_monotonic(self, enriched_dataset):
        """Each step's cumulative_yield is <= the previous step."""
        planner = SynthesisPlanner(enriched_dataset)
        plan = planner.plan_synthesis("A_AAA", modifications={})
        for i in range(1, len(plan.steps)):
            assert plan.steps[i].cumulative_yield <= plan.steps[i - 1].cumulative_yield


class TestPlanWithModification:
    """Tests for modified-position synthesis planning."""

    def test_plan_with_modification(self, enriched_dataset):
        """Modified position uses the correct modified monomer code."""
        planner = SynthesisPlanner(enriched_dataset)
        # Position 0 of A_AAA is A; apply 2OMe there
        plan = planner.plan_synthesis("A_AAA", modifications={0: "2OMe"})

        # Find the step for position 0 (it will be the last step since 3'->5')
        pos0_step = None
        for step in plan.steps:
            if step.position == 0:
                pos0_step = step
                break

        assert pos0_step is not None
        assert pos0_step.modification == "2OMe"
        assert "2OMe" in pos0_step.monomer_name

    def test_plan_incompatibility(self, enriched_dataset):
        """LNA + PS on adjacent positions produces an incompatibility warning."""
        planner = SynthesisPlanner(enriched_dataset)
        # Apply LNA at position 10, PS at position 11
        plan = planner.plan_synthesis("A_AAA", modifications={10: "LNA", 11: "PS"})
        assert len(plan.incompatibilities) > 0
        # The warning should mention the adjacent modifications
        warning_text = " ".join(plan.incompatibilities)
        assert "LNA" in warning_text
        assert "PS" in warning_text
