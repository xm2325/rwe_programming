from __future__ import annotations

import json
from pathlib import Path

from .deliverables import write_deliverables
from .independent_cox import reconcile_cox
from .longitudinal import build_analysis_cohort_python, make_longitudinal_source_population
from .longitudinal_qc import write_longitudinal_qc_manifest
from .longitudinal_validation import reconcile_longitudinal_builders
from .metadata import data_dictionary_frame
from .qc import cohort_attrition, write_qc_manifest
from .report import build_report
from .source_population import make_synthetic_source_population
from .study_config import DEFAULT_CONFIG, StudyConfig


def build_reviewer_bundle(
    output_dir: str | Path = "artifacts/reviewer_bundle",
    config: StudyConfig = DEFAULT_CONFIG,
    n_boot: int = 50,
) -> dict[str, str]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    paths: dict[str, str] = {}
    paths.update(write_deliverables(out, n=config.n_patients, seed=config.seed))
    paths["study_config"] = str(config.write_json(out / "study_config.json"))

    source = make_synthetic_source_population(n_eligible=config.n_patients, seed=config.seed)
    attrition_path = out / "cohort_attrition.csv"
    cohort_attrition(source, config).to_csv(attrition_path, index=False)
    paths["cohort_attrition"] = str(attrition_path)

    dictionary_path = out / "analysis_data_dictionary.csv"
    data_dictionary_frame().to_csv(dictionary_path, index=False)
    paths["analysis_data_dictionary"] = str(dictionary_path)

    paths["qc_manifest"] = str(write_qc_manifest(out / "qc_manifest.json", config))

    reconciliation_path = out / "independent_cox_reconciliation.json"
    reconciliation_path.write_text(
        json.dumps(reconcile_cox(n=config.n_patients, seed=config.seed), indent=2),
        encoding="utf-8",
    )
    paths["independent_cox_reconciliation"] = str(reconciliation_path)

    longitudinal_dir = out / "longitudinal_sources"
    longitudinal_dir.mkdir(exist_ok=True)
    sources = make_longitudinal_source_population(n_eligible=config.n_patients, seed=config.seed)
    source_inventory = {}
    for name, frame in sources.items():
        source_path = longitudinal_dir / f"{name}.csv"
        frame.to_csv(source_path, index=False)
        paths[f"source_{name}"] = str(source_path)
        source_inventory[name] = {"rows": int(len(frame)), "columns": list(frame.columns)}

    longitudinal_cohort = build_analysis_cohort_python(sources)
    longitudinal_cohort_path = out / "longitudinal_analysis_cohort.csv"
    longitudinal_cohort.to_csv(longitudinal_cohort_path, index=False)
    paths["longitudinal_analysis_cohort"] = str(longitudinal_cohort_path)

    source_inventory_path = out / "longitudinal_source_inventory.json"
    source_inventory_path.write_text(json.dumps(source_inventory, indent=2), encoding="utf-8")
    paths["longitudinal_source_inventory"] = str(source_inventory_path)

    paths["longitudinal_qc_manifest"] = str(
        write_longitudinal_qc_manifest(
            out / "longitudinal_qc_manifest.json",
            sources=sources,
            n=config.n_patients,
            seed=config.seed,
        )
    )

    longitudinal_reconciliation_path = out / "longitudinal_builder_reconciliation.json"
    longitudinal_reconciliation_path.write_text(
        json.dumps(
            reconcile_longitudinal_builders(
                n=config.n_patients,
                seed=config.seed,
                sql_path="sql/longitudinal_analysis_cohort.sql",
            ),
            indent=2,
        ),
        encoding="utf-8",
    )
    paths["longitudinal_builder_reconciliation"] = str(longitudinal_reconciliation_path)

    report_path = out / "rwe_validation_report.html"
    build_report(report_path, n=config.n_patients, seed=config.seed, n_boot=n_boot)
    paths["validation_report"] = str(report_path)

    inventory_path = out / "bundle_inventory.json"
    paths["bundle_inventory"] = str(inventory_path)
    inventory_path.write_text(json.dumps(paths, indent=2, sort_keys=True), encoding="utf-8")
    return paths


def main() -> None:
    print(json.dumps(build_reviewer_bundle(), indent=2))


if __name__ == "__main__":
    main()
