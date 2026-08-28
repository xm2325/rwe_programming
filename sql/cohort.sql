-- Illustrative cohort construction for the restored synthetic RWE portfolio.
-- This is intentionally schema-neutral and uses no proprietary data.
WITH eligible AS (
    SELECT *
    FROM person_level_source
    WHERE age >= 18
      AND baseline_egfr >= 45
      AND prevalent_ckd = 0
),
index_cohort AS (
    SELECT
        patient_id,
        index_date,
        CASE WHEN uncontrolled_gout = 1 THEN 1 ELSE 0 END AS ucg,
        age,
        female,
        diabetes,
        hypertension,
        baseline_egfr AS egfr,
        baseline_urate,
        prior_flares
    FROM eligible
)
SELECT * FROM index_cohort;
