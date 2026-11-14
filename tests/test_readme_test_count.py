"""The README's test count must equal the tests that actually exist.

The number in "Layout" is hand-typed, so it goes stale the moment a test is
added, and a reader uses it to judge whether the suite is substantial.
Counting `def test_` is not enough because parametrised tests expand, so this
counts what pytest actually collects, in a subprocess.
"""
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
README = (ROOT / "README.md").read_text()


def _collected() -> int:
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", "tests", "-q", "--collect-only"],
        cwd=ROOT, capture_output=True, text=True, timeout=600)
    out = proc.stdout + proc.stderr
    m = re.search(r"(\d+) tests? collected", out)
    if not m:
        pytest.fail(f"could not read a collected count from pytest:\n{out[-2000:]}")
    return int(m.group(1))


def test_readme_test_count_matches_the_collected_count():
    """`tests/   N tests` in the Layout section must be the real N."""
    m = re.search(r"^tests/\s+(\d+)\s+tests", README, re.M)
    assert m, "the Layout section no longer states a test count in the " \
              "form 'tests/   N tests'; update this test if the wording moved"
    stated = int(m.group(1))
    actual = _collected()
    assert stated == actual, (
        f"README says {stated} tests, pytest collects {actual}. "
        f"Update the README count.")
