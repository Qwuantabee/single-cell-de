"""The claim auditor must pass on the README and must fail when a number in it moves.

A checker that never fails is not a check, so the second test edits one quoted
number in a copy of the README and requires the auditor to reject it.
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "scripts" / "audit_claims.py"


def _run(*args):
    return subprocess.run([sys.executable, str(AUDIT), "--strict", *args],
                          cwd=ROOT, capture_output=True, text=True, timeout=300)


def test_every_readme_number_matches_its_artefact():
    proc = _run()
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "all anchored to a phrase" in proc.stdout


def test_auditor_rejects_a_moved_number(tmp_path):
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    # the headline power point: pseudobulk found 48.8 of 50 true genes
    needle = "with 8 donors it found 48.8 of 50"
    assert needle in readme, "headline sentence moved; update this test's needle"
    mutated = tmp_path / "README.md"
    mutated.write_text(readme.replace(needle, needle.replace("48.8", "49.8")),
                       encoding="utf-8")
    proc = _run("--readme", str(mutated))
    assert proc.returncode == 1, proc.stdout
    assert "headline pseudobulk TP" in proc.stdout
