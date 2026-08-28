from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from rwe_programming.nhanes_ps import (
    fit_nhanes_propensity_score,
    nhanes_ps_diagnostics,
    prepare_nhanes_gout_ps,
)


def read_xpt(path: Path) -> pd.DataFrame:
    return pd.read_sas(path, format="xport", encoding="utf-8")


def main(
    data_dir: str = "data/nhanes_2017_2018",
    output_dir: str = "artifacts/nhanes_ps",
) -> None:
    data = Path(data_dir)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    frame = prepare_nhanes_gout_ps(
        demo=read_xpt(data / "DEMO_J.XPT"),
        mcq=read_xpt(data / "MCQ_J.XPT"),
        rx=read_xpt(data / "RXQ_RX_J.XPT"),
        biopro=read_xpt(data / "BIOPRO_J.XPT"),
        bmx=read_xpt(data / "BMX_J.XPT"),
        diabetes=read_xpt(data / "DIQ_J.XPT"),
        blood_pressure=read_xpt(data / "BPQ_J.XPT"),
    )
    weighted = fit_nhanes_propensity_score(frame)
    diagnostics = nhanes_ps_diagnostics(weighted)

    # Keep only analysis-ready public NHANES fields; no synthetic substitution is made.
    weighted.to_csv(out / "nhanes_gout_ps_analysis.csv", index=False)
    (out / "nhanes_ps_diagnostics.json").write_text(
        json.dumps(diagnostics, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(diagnostics, indent=2))


if __name__ == "__main__":
    main()
