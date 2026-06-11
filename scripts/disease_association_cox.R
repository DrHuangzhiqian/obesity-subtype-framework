suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(survival)
  library(broom)
  library(openxlsx)
})

config <- list(
  analysis_data = "data/analysis_covariates.csv",
  outcome_dir = "data/outcomes",
  output_dir = "results/disease_association",
  id_col = "participant_id",
  subtype_col = "subtype",
  reference_subtype = "reference",
  time_col = "time_to_event",
  status_col = "event_status",
  outcome_suffix = "_outcome.csv",
  outcomes = c(
    "myocardial_infarction",
    "heart_failure",
    "stroke",
    "diabetic_microvascular_complications",
    "chronic_kidney_disease",
    "chronic_liver_disease"
  ),
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
  ),
  p_adjust_method = "BH"
)

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

analysis_data <- read_csv(config$analysis_data, show_col_types = FALSE) %>%
  mutate(
    "{config$subtype_col}" := factor(
      .data[[config$subtype_col]],
      levels = c(config$reference_subtype, setdiff(unique(.data[[config$subtype_col]]), config$reference_subtype))
    )
  )

fit_one_outcome <- function(outcome_name) {
  outcome_path <- file.path(config$outcome_dir, paste0(outcome_name, config$outcome_suffix))
  outcome_data <- read_csv(outcome_path, show_col_types = FALSE)

  required_cols <- c(config$id_col, config$time_col, config$status_col)
  missing_cols <- setdiff(required_cols, names(outcome_data))
  if (length(missing_cols) > 0) {
    stop("Missing required outcome columns for ", outcome_name, ": ", paste(missing_cols, collapse = ", "))
  }

  model_data <- analysis_data %>%
    inner_join(outcome_data, by = config$id_col) %>%
    filter(.data[[config$time_col]] > 0) %>%
    select(all_of(c(config$subtype_col, config$time_col, config$status_col, config$covariates))) %>%
    filter(complete.cases(.))

  surv_obj <- Surv(time = model_data[[config$time_col]], event = model_data[[config$status_col]])
  model_formula <- as.formula(
    paste("surv_obj ~", paste(c(config$subtype_col, config$covariates), collapse = " + "))
  )

  cox_fit <- coxph(model_formula, data = model_data)

  tidy(cox_fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(grepl(paste0("^", config$subtype_col), term)) %>%
    transmute(
      outcome = outcome_name,
      term,
      n = nrow(model_data),
      events = sum(model_data[[config$status_col]] == 1, na.rm = TRUE),
      HR = estimate,
      CI_low = conf.low,
      CI_high = conf.high,
      p_value = p.value
    )
}

cox_results <- map_dfr(config$outcomes, fit_one_outcome) %>%
  group_by(term) %>%
  mutate(p_adjusted = p.adjust(p_value, method = config$p_adjust_method)) %>%
  ungroup() %>%
  mutate(
    HR_95CI = sprintf("%.2f (%.2f, %.2f)", HR, CI_low, CI_high),
    p_value_display = if_else(p_value < 0.001, "<0.001", sprintf("%.3f", p_value)),
    p_adjusted_display = if_else(p_adjusted < 0.001, "<0.001", sprintf("%.3f", p_adjusted))
  )

write_csv(cox_results, file.path(config$output_dir, "cox_subtype_outcome_results.csv"))

workbook <- createWorkbook()
for (outcome_name in unique(cox_results$outcome)) {
  addWorksheet(workbook, sheetName = substr(outcome_name, 1, 31))
  writeData(workbook, sheet = substr(outcome_name, 1, 31), cox_results %>% filter(outcome == outcome_name))
}
saveWorkbook(workbook, file.path(config$output_dir, "cox_subtype_outcome_results.xlsx"), overwrite = TRUE)
