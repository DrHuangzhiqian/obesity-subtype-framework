suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(glmnet)
  library(caret)
  library(pROC)
})

config <- list(
  subtype_data = "data/subtype_assignments.csv",
  covariate_data = "data/covariates.csv",
  molecular_data = "data/molecular_features.csv",
  feature_screen = "data/significant_molecular_features.csv",
  output_dir = "results/molecular_signature",
  id_col = "participant_id",
  subtype_col = "subtype",
  reference_subtype = "reference",
  comparison_subtype = "comparison",
  feature_col = "feature",
  contrast_col = "contrast",
  adjusted_p_col = "adjusted_p",
  contrast_label = "comparison_vs_reference",
  adjusted_p_threshold = 0.05,
  train_fraction = 0.70,
  seed = 1234,
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
    "antihypertensive_medication",
    "fasting_time"
  )
)

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(config$seed)

subtype_data <- read_csv(config$subtype_data, show_col_types = FALSE) %>%
  select(all_of(c(config$id_col, config$subtype_col))) %>%
  filter(.data[[config$subtype_col]] %in% c(config$reference_subtype, config$comparison_subtype)) %>%
  mutate(
    outcome_binary = as.integer(.data[[config$subtype_col]] == config$comparison_subtype)
  )

covariate_data <- read_csv(config$covariate_data, show_col_types = FALSE)
molecular_data <- read_csv(config$molecular_data, show_col_types = FALSE)
feature_screen <- read_csv(config$feature_screen, show_col_types = FALSE)

candidate_features <- feature_screen %>%
  filter(
    .data[[config$contrast_col]] == config$contrast_label,
    .data[[config$adjusted_p_col]] < config$adjusted_p_threshold
  ) %>%
  pull(all_of(config$feature_col)) %>%
  unique() %>%
  intersect(names(molecular_data))

if (length(candidate_features) == 0) {
  stop("No candidate molecular features were available after screening.")
}

model_data <- subtype_data %>%
  inner_join(covariate_data, by = config$id_col) %>%
  inner_join(molecular_data %>% select(all_of(c(config$id_col, candidate_features))), by = config$id_col) %>%
  filter(complete.cases(select(., all_of(c("outcome_binary", config$covariates, candidate_features)))))

train_index <- createDataPartition(model_data$outcome_binary, p = config$train_fraction, list = FALSE)
train_data <- model_data[train_index, ]
validation_data <- model_data[-train_index, ]

formula_str <- paste("~", paste(c(config$covariates, candidate_features), collapse = " + "))
x_train <- model.matrix(as.formula(formula_str), train_data)[, -1, drop = FALSE]
y_train <- train_data$outcome_binary

x_covariates <- model.matrix(
  as.formula(paste("~", paste(config$covariates, collapse = " + "))),
  train_data
)[, -1, drop = FALSE]

penalty_factor <- ifelse(colnames(x_train) %in% colnames(x_covariates), 0, 1)

cv_fit <- cv.glmnet(
  x = x_train,
  y = y_train,
  family = "binomial",
  alpha = 1,
  penalty.factor = penalty_factor,
  nfolds = 10,
  type.measure = "deviance"
)

selected_lambda <- cv_fit$lambda.1se
coefficient_table <- as.matrix(coef(cv_fit, s = selected_lambda)) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("feature") %>%
  rename(coefficient = s1) %>%
  filter(feature != "(Intercept)", coefficient != 0)

molecular_coefficients <- coefficient_table %>%
  filter(feature %in% candidate_features)

write_csv(coefficient_table, file.path(config$output_dir, "lasso_coefficients_all.csv"))
write_csv(molecular_coefficients, file.path(config$output_dir, "lasso_coefficients_molecular.csv"))

score_data <- function(dat) {
  dat$signature_score <- as.numeric(
    as.matrix(dat[, molecular_coefficients$feature, drop = FALSE]) %*% molecular_coefficients$coefficient
  )
  dat
}

train_scored <- score_data(train_data) %>% mutate(dataset = "training")
validation_scored <- score_data(validation_data) %>% mutate(dataset = "validation")

signature_scores <- bind_rows(train_scored, validation_scored) %>%
  select(all_of(config$id_col), all_of(config$subtype_col), outcome_binary, dataset, signature_score)

write_csv(signature_scores, file.path(config$output_dir, "molecular_signature_scores.csv"))

evaluate_auc <- function(dat, label) {
  roc_obj <- roc(dat$outcome_binary, dat$signature_score, quiet = TRUE)
  auc_ci <- ci.auc(roc_obj)
  tibble::tibble(
    dataset = label,
    auc = as.numeric(auc(roc_obj)),
    auc_ci_low = as.numeric(auc_ci[1]),
    auc_ci_high = as.numeric(auc_ci[3])
  )
}

auc_results <- bind_rows(
  evaluate_auc(train_scored, "training"),
  evaluate_auc(validation_scored, "validation")
)

write_csv(auc_results, file.path(config$output_dir, "molecular_signature_auc.csv"))
