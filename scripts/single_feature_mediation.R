suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(mediation)
  library(openxlsx)
})

config <- list(
  analysis_data = "data/mediation_analysis_data.csv",
  mediation_plan = "data/mediation_plan.csv",
  output_dir = "results/mediation",
  subtype_col = "subtype",
  reference_subtype = "reference",
  comparison_subtype = "comparison",
  outcome_col = "outcome",
  mediator_col = "mediator",
  outcome_status_col = "event_status",
  sims = 1000,
  covariates = c(
    "age",
    "sex",
    "ethnicity",
    "socioeconomic_index",
    "education",
    "smoking_status",
    "alcohol_status",
    "lipid_lowering_medication",
    "antidiabetic_medication",
    "antihypertensive_medication"
  )
)

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

analysis_data <- read_csv(config$analysis_data, show_col_types = FALSE) %>%
  filter(.data[[config$subtype_col]] %in% c(config$reference_subtype, config$comparison_subtype)) %>%
  mutate(
    "{config$subtype_col}" := factor(
      .data[[config$subtype_col]],
      levels = c(config$reference_subtype, config$comparison_subtype)
    )
  )

mediation_plan <- read_csv(config$mediation_plan, show_col_types = FALSE)

required_plan_cols <- c(config$outcome_col, config$mediator_col)
missing_plan_cols <- setdiff(required_plan_cols, names(mediation_plan))
if (length(missing_plan_cols) > 0) {
  stop("Missing required mediation plan columns: ", paste(missing_plan_cols, collapse = ", "))
}

run_one_mediation <- function(outcome_name, mediator_name) {
  required_cols <- c(config$subtype_col, config$outcome_status_col, mediator_name, config$covariates)
  missing_cols <- setdiff(required_cols, names(analysis_data))
  if (length(missing_cols) > 0) {
    return(tibble::tibble(
      outcome = outcome_name,
      mediator = mediator_name,
      status = "failed",
      reason = paste("Missing columns:", paste(missing_cols, collapse = ", "))
    ))
  }

  model_data <- analysis_data %>%
    filter(.data[[config$outcome_col]] == outcome_name) %>%
    select(all_of(required_cols)) %>%
    filter(complete.cases(.))

  if (nrow(model_data) < 50 || length(unique(model_data[[config$outcome_status_col]])) < 2) {
    return(tibble::tibble(
      outcome = outcome_name,
      mediator = mediator_name,
      status = "failed",
      reason = "Insufficient complete cases or outcome variation"
    ))
  }

  mediator_formula <- as.formula(
    paste0("`", mediator_name, "` ~ ", config$subtype_col, " + ", paste(config$covariates, collapse = " + "))
  )
  outcome_formula <- as.formula(
    paste0(config$outcome_status_col, " ~ ", config$subtype_col, " + `", mediator_name, "` + ", paste(config$covariates, collapse = " + "))
  )

  mediator_model <- lm(mediator_formula, data = model_data)
  outcome_model <- glm(outcome_formula, family = binomial(link = "logit"), data = model_data)

  med_fit <- tryCatch(
    mediate(
      model.m = mediator_model,
      model.y = outcome_model,
      treat = config$subtype_col,
      mediator = mediator_name,
      boot = TRUE,
      sims = config$sims
    ),
    error = function(e) e
  )

  if (inherits(med_fit, "error")) {
    return(tibble::tibble(
      outcome = outcome_name,
      mediator = mediator_name,
      status = "failed",
      reason = conditionMessage(med_fit)
    ))
  }

  med_summary <- summary(med_fit)

  tibble::tibble(
    outcome = outcome_name,
    mediator = mediator_name,
    status = "ok",
    n = nrow(model_data),
    ACME = med_summary$d0,
    ACME_ci_low = med_summary$d0.ci[1],
    ACME_ci_high = med_summary$d0.ci[2],
    ACME_p = med_summary$d0.p,
    ADE = med_summary$z0,
    ADE_ci_low = med_summary$z0.ci[1],
    ADE_ci_high = med_summary$z0.ci[2],
    ADE_p = med_summary$z0.p,
    total_effect = med_summary$tau.coef,
    total_effect_ci_low = med_summary$tau.ci[1],
    total_effect_ci_high = med_summary$tau.ci[2],
    total_effect_p = med_summary$tau.p,
    proportion_mediated = med_summary$n0,
    proportion_mediated_ci_low = med_summary$n0.ci[1],
    proportion_mediated_ci_high = med_summary$n0.ci[2],
    proportion_mediated_p = med_summary$n0.p
  )
}

mediation_results <- pmap_dfr(
  list(
    mediation_plan[[config$outcome_col]],
    mediation_plan[[config$mediator_col]]
  ),
  run_one_mediation
)

write_csv(mediation_results, file.path(config$output_dir, "single_feature_mediation_results.csv"))
write.xlsx(
  list(
    mediation_results = mediation_results,
    mediation_plan = mediation_plan
  ),
  file.path(config$output_dir, "single_feature_mediation_results.xlsx"),
  overwrite = TRUE
)
