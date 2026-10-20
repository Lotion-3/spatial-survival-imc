"""Run the R tests (tests/R/): the Plevritis-lab CLQ code as driven by scripts/clq_official.R, and the
output format of the CELESTA bridge scripts/celesta_official.R. Skipped when R or the R packages are
not installed (CI has no R); run scripts/setup_r.R first."""

import glob
import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _rscript() -> str | None:
    found = shutil.which("Rscript")
    if found:
        return found
    hits = sorted(glob.glob(r"C:\Program Files\R\R-*\bin\Rscript.exe"))
    return hits[-1] if hits else None


def test_r_suite():
    rscript = _rscript()
    if rscript is None:
        pytest.skip("Rscript not found")
    if not (ROOT / "third_party" / "Spatial_Permutation_and_Normalization").exists():
        pytest.skip("run scripts/setup_r.R first (clones the lab repositories)")
    probe = subprocess.run([rscript, "-e", "library(CELESTA); library(testthat)"], capture_output=True)
    if probe.returncode != 0:
        pytest.skip("CELESTA / testthat not installed; run scripts/setup_r.R")
    r = subprocess.run([rscript, str(ROOT / "tests" / "R" / "run_tests.R")], capture_output=True, text=True,
                       env=os.environ | {"SPATIALSURV_ROOT": str(ROOT)}, timeout=1800)
    assert r.returncode == 0, r.stdout[-4000:] + r.stderr[-4000:]
