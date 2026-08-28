from __future__ import annotations

from pathlib import Path
import html
import json

from .pipeline import run_pipeline
from .sensitivity import bootstrap_hr, proportional_hazards_diagnostic
from .validation import (
    missingness_sensitivity,
    negative_control_analysis,
    omop_shape_reconciliation,
    outcome_sensitivity,
    sql_pandas_reconciliation,
    weight_trimming_sensitivity,
)


def build_report(
    path: str | Path = "validation/rwe_validation_report.html",
    n_boot: int = 30,
) -> Path:
    sections = {
        "primary": run_pipeline(),
        "sql_pandas": sql_pandas_reconciliation(),
        "omop_shape": omop_shape_reconciliation(),
        "weight_trimming": weight_trimming_sensitivity(),
        "missingness": missingness_sensitivity(),
        "bootstrap": bootstrap_hr(n_boot=n_boot),
        "ph_diagnostic": proportional_hazards_diagnostic(),
        "negative_control": negative_control_analysis(),
        "outcome_sensitivity": outcome_sensitivity(),
    }
    rows = "".join(
        f"<h2>{html.escape(name)}</h2><pre>{html.escape(json.dumps(result, indent=2))}</pre>"
        for name, result in sections.items()
    )
    document = f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>RWE validation report</title>
<style>
body {{ font-family: system-ui, sans-serif; max-width: 1000px; margin: 40px auto; line-height: 1.45; }}
pre {{ background: #f5f5f5; padding: 14px; overflow: auto; }}
</style>
</head>
<body>
<h1>Auditable RWE validation report</h1>
<p>All patient-level data in this report are synthetic.</p>
{rows}
</body>
</html>"""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(document, encoding="utf-8")
    return output
