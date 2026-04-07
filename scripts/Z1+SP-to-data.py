#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from z1plus_tools.cli import main_sp_to_data


if __name__ == '__main__':
    raise SystemExit(main_sp_to_data())
