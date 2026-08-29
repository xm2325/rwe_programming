from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import json


@dataclass(frozen=True)
class StudyConfig:
    study_id: str = "RWE-GOUT-CKD-001"
    study_version: str = "0.10.0-restored"
    population: str = "synthetic adults with gout"
    exposure: str = "uncontrolled gout"
    comparator: str = "controlled gout"
    primary_outcome: str = "incident CKD"
    min_age: int = 18
    min_baseline_egfr: float = 45.0
    max_followup_years: float = 5.0
    propensity_clip_low: float = 0.01
    propensity_clip_high: float = 0.99
    weight_trim_low_quantile: float = 0.01
    weight_trim_high_quantile: float = 0.99
    balance_threshold_abs_smd: float = 0.10
    seed: int = 20260817
    n_patients: int = 9184

    def to_dict(self) -> dict:
        return asdict(self)

    def write_json(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")
        return path


DEFAULT_CONFIG = StudyConfig()
