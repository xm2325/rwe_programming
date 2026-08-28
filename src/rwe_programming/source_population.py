from __future__ import annotations

import numpy as np
import pandas as pd

from .pipeline import make_synthetic_cohort
from .study_config import DEFAULT_CONFIG, StudyConfig


def make_synthetic_source_population(
    n_eligible: int = DEFAULT_CONFIG.n_patients,
    seed: int = DEFAULT_CONFIG.seed,
    ineligible_fraction: float = 0.12,
) -> pd.DataFrame:
    """Create a synthetic source population with deterministic staged exclusions.

    The requested ``n_eligible`` remains the final analysis-cohort size. Extra records
    are added solely to exercise source-to-analysis eligibility logic and QC.
    """
    eligible = make_synthetic_cohort(n=n_eligible, seed=seed).copy()
    n_extra = max(4, int(round(n_eligible * ineligible_fraction)))
    extras = make_synthetic_cohort(n=n_extra, seed=seed + 1).copy()
    extras["patient_id"] = np.arange(n_eligible + 1, n_eligible + n_extra + 1)

    groups = np.array_split(np.arange(n_extra), 4)
    # Ensure each record violates exactly the intended first applicable rule.
    if len(groups[0]):
        extras.loc[groups[0], "age"] = 17.0
    if len(groups[1]):
        extras.loc[groups[1], "age"] = np.maximum(extras.loc[groups[1], "age"], 18.0)
        extras.loc[groups[1], "egfr"] = 30.0
    if len(groups[2]):
        extras.loc[groups[2], "age"] = np.maximum(extras.loc[groups[2], "age"], 18.0)
        extras.loc[groups[2], "egfr"] = np.maximum(extras.loc[groups[2], "egfr"], 45.0)
        extras.loc[groups[2], "followup_years"] = 0.0
    if len(groups[3]):
        extras.loc[groups[3], "age"] = np.maximum(extras.loc[groups[3], "age"], 18.0)
        extras.loc[groups[3], "egfr"] = np.maximum(extras.loc[groups[3], "egfr"], 45.0)
        extras.loc[groups[3], "followup_years"] = np.maximum(extras.loc[groups[3], "followup_years"], 1e-6)
        extras.loc[groups[3], "ucg"] = 2

    return pd.concat([eligible, extras], ignore_index=True)


def select_analysis_cohort(
    source: pd.DataFrame,
    config: StudyConfig = DEFAULT_CONFIG,
) -> pd.DataFrame:
    out = source.copy()
    out = out[out["age"] >= config.min_age]
    out = out[out["egfr"] >= config.min_baseline_egfr]
    out = out[out["followup_years"] > 0]
    out = out[out["ucg"].isin([0, 1])]
    return out.sort_values("patient_id").reset_index(drop=True)
