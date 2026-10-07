"""Model definitions. All fitted preprocessing (imputation, scaling) lives inside the
pipeline so that it is re-fit on each training fold."""

from __future__ import annotations

import warnings
from contextlib import contextmanager
from dataclasses import dataclass, field

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.compose import ColumnTransformer
from sklearn.exceptions import ConvergenceWarning
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline, make_pipeline
from sklearn.preprocessing import StandardScaler
from sksurv.linear_model import CoxnetSurvivalAnalysis, CoxPHSurvivalAnalysis

from .features import clr


class CLRTransformer(BaseEstimator, TransformerMixin):
    """Row-wise CLR of count columns. Stateless: fit() learns nothing from the data."""

    def __init__(self, pseudocount: float = 0.5):
        self.pseudocount = pseudocount

    def fit(self, X, y=None):
        self.n_features_in_ = np.asarray(X).shape[1]
        return self

    def transform(self, X):
        return clr(np.asarray(X, dtype=float), self.pseudocount)


class ResidualizeSpatial(BaseEstimator, TransformerMixin):
    """Replace spatial features by their residuals from a linear regression on
    [CLR(composition counts), log(total cells)], fit on training rows only.

    Input columns: the n_spatial spatial features first, then the composition counts.
    Spatial NaNs are median-imputed (training medians) before the regression.
    Isolates "arrangement beyond amount": the part of each spatial feature that
    composition and image size do not explain linearly.
    """

    def __init__(self, n_spatial: int, pseudocount: float = 0.5):
        self.n_spatial = n_spatial
        self.pseudocount = pseudocount

    def _design(self, X):
        X = np.asarray(X, dtype=float)
        S, C = X[:, : self.n_spatial], X[:, self.n_spatial :]
        Z = np.c_[clr(C, self.pseudocount), np.log(C.sum(axis=1) + 1.0)]
        return S, Z

    def fit(self, X, y=None):
        S, Z = self._design(X)
        self.medians_ = np.nanmedian(S, axis=0)
        self.medians_ = np.where(np.isnan(self.medians_), 0.0, self.medians_)
        S = np.where(np.isnan(S), self.medians_, S)
        self.s_mean_, self.z_mean_ = S.mean(axis=0), Z.mean(axis=0)
        # CLR columns sum to zero (rank-deficient); lstsq returns the minimum-norm solution.
        self.coef_ = np.linalg.lstsq(Z - self.z_mean_, S - self.s_mean_, rcond=None)[0]
        return self

    def transform(self, X):
        S, Z = self._design(X)
        S = np.where(np.isnan(S), self.medians_, S)
        return S - self.s_mean_ - (Z - self.z_mean_) @ self.coef_


@dataclass
class FeatureSet:
    name: str
    clinical: list[str]
    composition: list[str] = field(default_factory=list)
    spatial: list[str] = field(default_factory=list)
    residualize_spatial: bool = False  # part II section B

    @property
    def columns(self) -> list[str]:
        return self.clinical + self.composition + self.spatial


def make_preprocessor(fs: FeatureSet, pseudocount: float) -> ColumnTransformer:
    blocks = [("clin", make_pipeline(SimpleImputer(strategy="median", keep_empty_features=True), StandardScaler()), fs.clinical)]
    if fs.composition:
        blocks.append(("comp", make_pipeline(CLRTransformer(pseudocount), StandardScaler()), fs.composition))
    if fs.spatial and fs.residualize_spatial:
        if not fs.composition:
            raise ValueError("residualising spatial features needs composition columns")
        blocks.append(("spat", make_pipeline(ResidualizeSpatial(len(fs.spatial), pseudocount), StandardScaler()),
                       fs.spatial + fs.composition))
    elif fs.spatial:
        blocks.append(("spat", make_pipeline(SimpleImputer(strategy="median", keep_empty_features=True), StandardScaler()), fs.spatial))
    # Column order of the output = clinical, composition, spatial (relied on for penalty_factor).
    return ColumnTransformer(blocks, remainder="drop", verbose_feature_names_out=False)


def make_model(fs: FeatureSet, cfg: dict, alpha: float | np.ndarray | None = None) -> Pipeline:
    """Clinical-only: ridge Cox (CoxPHSurvivalAnalysis, alpha tuned).
    With omics blocks: elastic-net Coxnet, clinical columns unpenalized (penalty_factor 0)."""
    pre = make_preprocessor(fs, cfg["clr_pseudocount"])
    if not fs.composition and not fs.spatial:
        est = CoxPHSurvivalAnalysis(alpha=0.0 if alpha is None else float(alpha), ties="breslow", n_iter=200)
    else:
        pf = np.r_[np.zeros(len(fs.clinical)), np.ones(len(fs.composition) + len(fs.spatial))]
        kw = dict(
            l1_ratio=cfg["l1_ratio"],
            penalty_factor=pf,
            max_iter=cfg["coxnet_max_iter"],
            normalize=False,
            fit_baseline_model=False,
        )
        if alpha is None:
            est = CoxnetSurvivalAnalysis(n_alphas=cfg["n_alphas"], alpha_min_ratio=cfg["alpha_min_ratio"], **kw)
        else:
            est = CoxnetSurvivalAnalysis(alphas=np.atleast_1d(alpha), **kw)
    return Pipeline([("pre", pre), ("cox", est)])


def is_coxnet(fs: FeatureSet) -> bool:
    return bool(fs.composition or fs.spatial)


_FIT_ERRORS = (ValueError, np.linalg.LinAlgError, ArithmeticError)


@contextmanager
def _quiet():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", (ConvergenceWarning, UserWarning, RuntimeWarning))
        yield


def _fit_coxnet_path(fs: FeatureSet, X, y, alphas: np.ndarray, cfg: dict) -> tuple[Pipeline, int]:
    """Fit Coxnet along a decreasing alpha path. If coxnet diverges (it can for the least
    penalised end of the path when events are few), drop the smallest alphas until it
    fits. Returns (model, number of leading alphas that were fitted)."""
    alphas = np.sort(np.asarray(alphas))[::-1]
    for n_ok in range(len(alphas), 0, -1):
        m = make_model(fs, cfg, alpha=alphas[:n_ok])
        try:
            m.fit(X, y)
            return m, n_ok
        except ArithmeticError:  # coxnet divergence only; other errors are bugs and must surface
            continue
    raise RuntimeError("Coxnet failed to fit even at the largest alpha")


def alpha_grid(fs: FeatureSet, X, y, cfg: dict) -> np.ndarray:
    """Penalty grid, derived from training data only. Same number of points for every model.

    Coxnet: the automatic path (n_alphas points from alpha_max down to
    alpha_min_ratio * alpha_max). If the full path diverges, alpha_min_ratio is raised
    (x2 each try) so the grid keeps n_alphas points. Ridge: fixed log-spaced grid."""
    if not is_coxnet(fs):
        return np.logspace(cfg["ridge_log10_min"], cfg["ridge_log10_max"], cfg["n_alphas"])[::-1]
    ratio = cfg["alpha_min_ratio"]
    with _quiet():
        while ratio < 1:
            m = make_model(fs, {**cfg, "alpha_min_ratio": ratio})
            try:
                m.fit(X, y)
                return np.asarray(m.named_steps["cox"].alphas_)
            except ArithmeticError:
                ratio *= 2
    raise RuntimeError("could not derive an alpha grid")


def fit_path(fs: FeatureSet, X, y, alphas: np.ndarray, cfg: dict) -> list:
    """Fit the model for every alpha. Returns one callable X -> risk per alpha, or None
    where the fit failed (those alphas are excluded from selection for this fold)."""
    with _quiet():
        if is_coxnet(fs):
            m, n_ok = _fit_coxnet_path(fs, X, y, alphas, cfg)
            ok = set(np.sort(np.asarray(alphas))[::-1][:n_ok])
            return [(lambda Xt, a=a, m=m: m.predict(Xt, alpha=a)) if a in ok else None for a in alphas]
        out = []
        for a in alphas:
            m = make_model(fs, cfg, alpha=a)
            try:
                m.fit(X, y)
                out.append(m.predict)
            except _FIT_ERRORS:
                out.append(None)
        return out


def fit_final(fs: FeatureSet, X, y, alpha: float, cfg: dict) -> Pipeline:
    """Refit on the full (outer) training set at the selected alpha. Coxnet is fitted along
    the training-data path down to alpha (warm starts); if that diverges, the smallest
    alpha that fits is used and recorded in chosen_alpha_."""
    with _quiet():
        if is_coxnet(fs):
            grid = alpha_grid(fs, X, y, cfg)
            grid = np.r_[grid[grid > alpha], alpha]
            m, n_ok = _fit_coxnet_path(fs, X, y, grid, cfg)
            m.chosen_alpha_ = float(grid[n_ok - 1])
            return m
        m = make_model(fs, cfg, alpha=alpha)
        m.fit(X, y)
        m.chosen_alpha_ = alpha
        return m


def predict_risk(model: Pipeline, X) -> np.ndarray:
    if isinstance(model.named_steps["cox"], CoxnetSurvivalAnalysis):
        return model.predict(X, alpha=model.chosen_alpha_)
    return model.predict(X)


def coefficients(model: Pipeline) -> np.ndarray:
    cox = model.named_steps["cox"]
    if isinstance(cox, CoxnetSurvivalAnalysis):
        j = int(np.argmin(np.abs(cox.alphas_ - model.chosen_alpha_)))
        return np.asarray(cox.coef_[:, j])
    return np.asarray(cox.coef_)
