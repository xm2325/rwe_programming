WITH base AS (
  SELECT p.patient_id,
         CAST(strftime('%Y', e.index_date) AS INTEGER) - p.birth_year AS age,
         p.female,
         e.index_date,
         e.enrollment_end,
         date(e.index_date, '-365 day') AS baseline_start,
         CASE
           WHEN e.enrollment_end < date(e.index_date, '+1825 day') THEN e.enrollment_end
           ELSE date(e.index_date, '+1825 day')
         END AS followup_end
  FROM patients p
  JOIN enrollment e USING (patient_id)
),
last_egfr AS (
  SELECT l.patient_id, l.value AS egfr
  FROM labs l
  JOIN base b USING (patient_id)
  WHERE l.lab_name = 'EGFR'
    AND l.lab_date >= b.baseline_start
    AND l.lab_date < b.index_date
    AND l.lab_date = (
      SELECT MAX(l2.lab_date) FROM labs l2
      WHERE l2.patient_id = l.patient_id
        AND l2.lab_name = 'EGFR'
        AND l2.lab_date >= b.baseline_start
        AND l2.lab_date < b.index_date
    )
),
last_urate AS (
  SELECT l.patient_id, l.value AS baseline_urate
  FROM labs l
  JOIN base b USING (patient_id)
  WHERE l.lab_name = 'SERUM_URATE'
    AND l.lab_date >= b.baseline_start
    AND l.lab_date < b.index_date
    AND l.lab_date = (
      SELECT MAX(l2.lab_date) FROM labs l2
      WHERE l2.patient_id = l.patient_id
        AND l2.lab_name = 'SERUM_URATE'
        AND l2.lab_date >= b.baseline_start
        AND l2.lab_date < b.index_date
    )
),
dx AS (
  SELECT b.patient_id,
         SUM(CASE WHEN d.code='GOUT_FLARE' THEN 1 ELSE 0 END) AS prior_flares,
         MAX(CASE WHEN d.code='DIABETES' THEN 1 ELSE 0 END) AS diabetes,
         MAX(CASE WHEN d.code='HYPERTENSION' THEN 1 ELSE 0 END) AS hypertension,
         MAX(CASE WHEN d.code='BASELINE_CKD' THEN 1 ELSE 0 END) AS baseline_ckd
  FROM base b
  LEFT JOIN diagnoses d ON d.patient_id=b.patient_id
    AND d.diagnosis_date >= b.baseline_start
    AND d.diagnosis_date < b.index_date
  GROUP BY b.patient_id
),
eligible AS (
  SELECT b.patient_id, b.age, b.female, b.index_date, b.followup_end,
         dx.diabetes, dx.hypertension, e.egfr, u.baseline_urate,
         dx.prior_flares,
         CASE WHEN u.baseline_urate >= 8.0 OR dx.prior_flares >= 2 THEN 1 ELSE 0 END AS ucg
  FROM base b
  JOIN last_egfr e USING (patient_id)
  JOIN last_urate u USING (patient_id)
  JOIN dx USING (patient_id)
  WHERE b.age >= 18
    AND e.egfr >= 45.0
    AND dx.baseline_ckd = 0
),
first_primary AS (
  SELECT o.patient_id, MIN(o.outcome_date) AS event_date
  FROM outcomes o
  JOIN eligible e USING (patient_id)
  WHERE o.outcome='INCIDENT_CKD'
    AND o.outcome_date > e.index_date
    AND o.outcome_date <= e.followup_end
  GROUP BY o.patient_id
),
first_strict AS (
  SELECT o.patient_id, MIN(o.outcome_date) AS event_date
  FROM outcomes o
  JOIN eligible e USING (patient_id)
  WHERE o.outcome='INCIDENT_CKD_CONFIRM'
    AND o.outcome_date > e.index_date
    AND o.outcome_date <= e.followup_end
  GROUP BY o.patient_id
),
first_negative AS (
  SELECT o.patient_id, MIN(o.outcome_date) AS event_date
  FROM outcomes o
  JOIN eligible e USING (patient_id)
  WHERE o.outcome='NEGATIVE_CONTROL'
    AND o.outcome_date > e.index_date
    AND o.outcome_date <= e.followup_end
  GROUP BY o.patient_id
)
SELECT e.patient_id, e.age, e.female, e.diabetes, e.hypertension,
       e.egfr, e.baseline_urate, e.prior_flares, e.ucg,
       (julianday(COALESCE(p.event_date, e.followup_end)) - julianday(e.index_date)) / 365.25 AS followup_years,
       CASE WHEN p.event_date IS NULL THEN 0 ELSE 1 END AS ckd_event,
       (julianday(COALESCE(s.event_date, e.followup_end)) - julianday(e.index_date)) / 365.25 AS ckd_strict_followup_years,
       CASE WHEN s.event_date IS NULL THEN 0 ELSE 1 END AS ckd_strict_event,
       (julianday(COALESCE(n.event_date, e.followup_end)) - julianday(e.index_date)) / 365.25 AS negative_followup_years,
       CASE WHEN n.event_date IS NULL THEN 0 ELSE 1 END AS negative_event
FROM eligible e
LEFT JOIN first_primary p USING (patient_id)
LEFT JOIN first_strict s USING (patient_id)
LEFT JOIN first_negative n USING (patient_id)
ORDER BY e.patient_id;
