from __future__ import annotations

from pathlib import Path
import html
import json

from .analysis import prepare_weighted_analysis, summarise_analysis
from .longitudinal import build_analysis_cohort_python, make_longitudinal_source_population
from .longitudinal_qc import longitudinal_qc_manifest
from .longitudinal_validation import reconcile_longitudinal_builders
from .pipeline import make_synthetic_cohort, run_pipeline
from .qc import cohort_attrition, qc_manifest
from .sensitivity import bootstrap_hr, proportional_hazards_diagnostic
from .source_population import make_synthetic_source_population
from .study_config import DEFAULT_CONFIG, StudyConfig
from .survival import survival_at_times, weighted_kaplan_meier
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
    n: int | None = None,
    seed: int | None = None,
) -> Path:
    config = StudyConfig(
        n_patients=DEFAULT_CONFIG.n_patients if n is None else n,
        seed=DEFAULT_CONFIG.seed if seed is None else seed,
    )
    raw = make_synthetic_cohort(n=config.n_patients, seed=config.seed)
    source = make_synthetic_source_population(n_eligible=config.n_patients, seed=config.seed)
    longitudinal_sources = make_longitudinal_source_population(n_eligible=config.n_patients, seed=config.seed)
    longitudinal_cohort = build_analysis_cohort_python(longitudinal_sources)
    weighted_longitudinal = prepare_weighted_analysis(longitudinal_cohort)
    km_summary = survival_at_times(weighted_kaplan_meier(weighted_longitudinal)).to_dict(orient="records")
    attrition = cohort_attrition(source, config).to_dict(orient="records")
    longitudinal_qc = longitudinal_qc_manifest(
        sources=longitudinal_sources,
        n=config.n_patients,
        seed=config.seed,
    )
    sections = {
        "study_configuration": config.to_dict(),
        "cohort_attrition": attrition,
        "runtime_qc_manifest": qc_manifest(config),
        "longitudinal_qc_manifest": longitudinal_qc,
        "primary_source_derived_analysis": summarise_analysis(longitudinal_cohort),
        "weighted_survival_summary": km_summary,
        "longitudinal_sql_python_reconciliation": reconcile_longitudinal_builders(
            n=config.n_patients,
            seed=config.seed,
            sql_path="sql/longitudinal_analysis_cohort.sql",
        ),
        "legacy_flat_reference_analysis": run_pipeline(n=config.n_patients, seed=config.seed),
        "flat_sql_pandas": sql_pandas_reconciliation(raw),
        "omop_shape": omop_shape_reconciliation(raw),
        "weight_trimming": weight_trimming_sensitivity(raw),
        "missingness": missingness_sensitivity(raw, seed=config.seed + 101),
        "bootstrap": bootstrap_hr(n_boot=n_boot, seed=config.seed + 202, df=raw),
        "ph_diagnostic": proportional_hazards_diagnostic(raw),
        "negative_control": negative_control_analysis(raw),
        "outcome_sensitivity": outcome_sensitivity(raw),
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
.pass {{ font-weight: 700; }}
</style>
</head>
<body>
<h1>Auditable RWE validation report</h1>
<p>Study: <strong>{html.escape(config.study_id)}</strong> · version {html.escape(config.study_version)}</p>
<p>All patient-level data in this report are synthetic. Runtime QC status: <span class="pass">{sections['runtime_qc_manifest']['status']}</span>; longitudinal QC status: <span class="pass">{sections['longitudinal_qc_manifest']['status']}</span>.</p>
<p>The primary analysis and weighted survival summaries are generated from the cohort derived from longitudinal source tables. The legacy flat cohort is retained only as a regression/reference validation layer.</p>
{rows}
</body>
</html>"""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(document, encoding="utf-8")
    return output
