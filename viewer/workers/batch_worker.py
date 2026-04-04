"""QThread wrapper that runs BatchScreener without blocking the Qt event loop."""

from __future__ import annotations

import logging

from PySide6.QtCore import QThread, Signal

from viewer.chemistry.batch_screener import BatchJob, BatchResult, BatchScreener

_log = logging.getLogger(__name__)


class BatchWorker(QThread):
    """Executes batch screening in a background thread.

    Qt delivers cross-thread signals via QueuedConnection by default, so
    result_ready always arrives on the main thread — safe to update widgets.
    """

    result_ready = Signal(object)        # BatchResult — fires per completed job
    progress_updated = Signal(int, int)  # (n_completed, n_total)
    batch_finished = Signal()

    def __init__(
        self,
        screener: BatchScreener,
        jobs: list[BatchJob],
        parent=None,
    ) -> None:
        super().__init__(parent)
        self._screener = screener
        self._jobs = jobs
        self._completed = 0

    def run(self) -> None:
        total = len(self._jobs)

        def on_result(result: BatchResult) -> None:
            self._completed += 1
            self.result_ready.emit(result)
            self.progress_updated.emit(self._completed, total)

        def on_done() -> None:
            self.batch_finished.emit()

        try:
            self._screener._run(self._jobs, on_result, on_done)
        except Exception:
            _log.exception("BatchWorker.run() failed")
            self.batch_finished.emit()

    def cancel(self) -> None:
        self._screener.cancel()
