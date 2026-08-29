args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) args[[1]] else "artifacts/cox_r_validation/weighted_analysis.csv"
output_path <- if (length(args) >= 2) args[[2]] else "artifacts/cox_r_validation/r_cox_reference.csv"

suppressPackageStartupMessages(library(survival))

d <- read.csv(input_path, stringsAsFactors = FALSE)
required <- c("patient_id", "followup_years", "ckd_event", "ucg", "stabilized_weight")
missing <- setdiff(required, names(d))
if (length(missing) > 0) {
  stop(paste("Missing required columns:", paste(missing, collapse = ", ")))
}

fit <- coxph(
  Surv(followup_years, ckd_event) ~ ucg,
  data = d,
  weights = stabilized_weight,
  ties = "breslow",
  robust = TRUE,
  cluster = patient_id,
  model = FALSE,
  x = FALSE,
  y = FALSE
)

out <- data.frame(
  implementation = "R survival::coxph",
  coefficient = as.numeric(coef(fit)[[1]]),
  hazard_ratio = exp(as.numeric(coef(fit)[[1]])),
  robust_se = sqrt(as.numeric(diag(vcov(fit))[[1]])),
  n = nrow(d),
  events = sum(d$ckd_event)
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.csv(out, output_path, row.names = FALSE)
print(out)
