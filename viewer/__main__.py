"""Entry point for the 3D pre-miRNA viewer application.

Run with: python -m viewer
"""

from __future__ import annotations

import logging
import sys


def main() -> None:
    logging.basicConfig(level=logging.WARNING, format="%(name)s: %(message)s")

    from PySide6.QtWidgets import QApplication

    from viewer.app import MainWindow

    app = QApplication.instance() or QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
