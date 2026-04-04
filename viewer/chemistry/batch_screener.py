"""Batch screening of multiple variants via a thread pool."""

from __future__ import annotations

import os
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Callable

if TYPE_CHECKING:
    from viewer.chemistry.double_screen import DoubleReplacementScreener, DoubleScreenResult
    from viewer.chemistry.virtual_screen import ScreenResult, VirtualScreener
    from viewer.data.schema import EnrichedVariantDataset


@dataclass
class BatchJob:
    """A single screening task for one variant."""

    variant_id: str
    mode: str       # "single" | "double"
    top_n: int = 20


@dataclass
class BatchResult:
    """Completed result for one BatchJob."""

    variant_id: str
    mode: str
    single_results: list = field(default_factory=list)   # list[ScreenResult]
    double_results: list = field(default_factory=list)   # list[DoubleScreenResult]
    error: str | None = None


class BatchScreener:
    """Runs screening jobs across multiple variants using a thread pool.

    Thread safety: screeners read immutable variant sequences and pre-fitted
    ridge weights; numpy operations release the GIL, so ThreadPoolExecutor
    gives real parallelism with no serialization overhead.
    """

    def __init__(
        self,
        screener: VirtualScreener,
        double_screener: DoubleReplacementScreener,
        dataset: EnrichedVariantDataset,
        worker_count: int | None = None,
    ) -> None:
        self._screener = screener
        self._double_screener = double_screener
        self._dataset = dataset
        self._worker_count = worker_count
        self._cancel = threading.Event()

    def cancel(self) -> None:
        self._cancel.set()

    def _run(
        self,
        jobs: list[BatchJob],
        on_result: Callable[[BatchResult], None],
        on_done: Callable[[], None],
    ) -> None:
        """Blocking execution; call from a background thread (e.g. BatchWorker).

        on_result is called from a worker thread — emit Qt signals there rather
        than touching widgets directly.
        """
        if not jobs:
            on_done()
            return

        self._cancel.clear()
        n_workers = self._worker_count or min(os.cpu_count() or 1, len(jobs))

        with ThreadPoolExecutor(max_workers=n_workers) as pool:
            futures = {pool.submit(self._run_one, job): job for job in jobs}
            for fut in as_completed(futures):
                if self._cancel.is_set():
                    break
                try:
                    result = fut.result()
                except Exception as exc:
                    job = futures[fut]
                    result = BatchResult(job.variant_id, job.mode, error=str(exc))
                on_result(result)

        on_done()

    def _run_one(self, job: BatchJob) -> BatchResult:
        if self._cancel.is_set():
            return BatchResult(job.variant_id, job.mode, error="cancelled")
        try:
            if job.mode == "single":
                hits = self._screener.rank_by_dc_ratio_shift(
                    job.variant_id, top_n=job.top_n
                )
                return BatchResult(job.variant_id, job.mode, single_results=hits)
            else:
                hits = self._double_screener.screen_double(
                    job.variant_id, top_n=job.top_n
                )
                return BatchResult(job.variant_id, job.mode, double_results=hits)
        except Exception as exc:
            return BatchResult(job.variant_id, job.mode, error=str(exc))
