"""Make `src/` importable in a fresh clone without `pip install -e .`.

Several test modules do `from scpipe import ...` with no path setup of their
own. The full suite passed only because an earlier-collected module inserted
`src/` into sys.path; running one of those files alone failed on collection
with ModuleNotFoundError.
"""

import sys
from pathlib import Path

SRC = Path(__file__).resolve().parent / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
