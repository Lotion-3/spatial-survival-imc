"""CELESTA-Lite: a vectorised Python port of the CELESTA cell-typing engine (marker
probabilities + spatial Potts/MRF refinement), used here to re-type the METABRIC IMC cells.

Independent re-implementation following the reference R code (plevritis-lab/CELESTA) and
Zhang et al., Nature Methods 2022. Not affiliated with the original authors.
"""

from __future__ import annotations

from .benchmark import compare_assignments
from .graph import SpatialGraph, delaunay_graph, knn_graph, neighbourhood_composition
from .mrf import UNKNOWN, CelestaConfig, CelestaResult, run_celesta
from .preprocessing import (
    MarkerModel,
    Signature,
    arcsinh_transform,
    artifact_mask,
    cell_type_scores,
    expression_probabilities,
    fit_marker_models,
)

__all__ = [
    "UNKNOWN",
    "CelestaConfig",
    "CelestaResult",
    "MarkerModel",
    "Signature",
    "SpatialGraph",
    "arcsinh_transform",
    "artifact_mask",
    "cell_type_scores",
    "compare_assignments",
    "delaunay_graph",
    "expression_probabilities",
    "fit_marker_models",
    "knn_graph",
    "neighbourhood_composition",
    "run_celesta",
]
