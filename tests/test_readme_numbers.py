"""Every markdown table row in README.md must appear verbatim in results/tables.md,
which scripts/report_tables.py generates from the saved results files."""

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_readme_tables_come_from_results():
    tables = ROOT / "results" / "tables.md"
    if not tables.exists():
        pytest.skip("run scripts/run.py and scripts/report_tables.py first")
    generated = set(tables.read_text(encoding="utf-8").splitlines())
    readme = (ROOT / "README.md").read_text(encoding="utf-8").splitlines()
    rows = [l for l in readme if l.startswith("| ") and not set(l) <= set("|- ")]
    # The data-files table in section 2 describes downloads, not results.
    rows = [r for r in rows if ".zip" not in r and not r.startswith("| File |")]
    # Layout tables (side-by-side figures) and text-only header rows carry no numbers.
    rows = [r for r in rows if "![" not in r and any(ch.isdigit() for ch in r)]
    assert rows, "no result tables found in README"
    missing = [r for r in rows if r not in generated]
    assert not missing, "README rows not in results/tables.md:\n" + "\n".join(missing)
