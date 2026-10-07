"""Evaluation helpers: marker-only vs. spatially refined assignments against ground truth."""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, classification_report, f1_score

__all__ = ["compare_assignments"]


def compare_assignments(truth: pd.Series, predictions: dict[str, pd.Series],
                        unassigned_label: str = "Unassigned") -> pd.DataFrame:
    """Accuracy, macro-F1 and coverage of several labelings against the same truth.

    Unassigned cells count as errors in ``accuracy`` and ``macro_f1``, and are excluded
    in ``accuracy_assigned``, so methods that abstain more are not rewarded for it.

    Args:
        truth: ground-truth label per cell.
        predictions: method name -> predicted label per cell (same index as ``truth``).
        unassigned_label: label used for cells left unassigned.

    Returns:
        One row per method.
    """
    rows = []
    for name, pred in predictions.items():
        if not pred.index.equals(truth.index):
            raise ValueError(f"prediction {name!r} is not aligned with the truth index")
        assigned = (pred != unassigned_label).to_numpy()
        labels = sorted(set(truth.unique()))
        rows.append(
            dict(
                method=name,
                coverage=float(assigned.mean()),
                accuracy=float(accuracy_score(truth, pred)),
                accuracy_assigned=float(accuracy_score(truth[assigned], pred[assigned])) if assigned.any() else np.nan,
                macro_f1=float(f1_score(truth, pred, labels=labels, average="macro", zero_division=0)),
            )
        )
    return pd.DataFrame(rows).set_index("method")


def report(truth: pd.Series, pred: pd.Series) -> str:
    """Per-class precision/recall/F1 (scikit-learn text report)."""
    return classification_report(truth, pred, zero_division=0)
