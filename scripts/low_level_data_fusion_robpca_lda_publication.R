# ============================================================
# LOW-LEVEL DATA FUSION WITH ROBPCA AND LDA
# ============================================================
# Purpose:
#   This script reproduces the low-level data fusion workflow used for
#   the classification of commercial cigarette samples based on inorganic,
#   organic, and isotopic descriptors.
#
# Main steps:
#   1. Data import and consistency checks
#   2. Low-level data fusion
#   3. ROBPCA-based multivariate outlier screening
#   4. Stratified train/test split
#   5. Block autoscaling and block scaling
#   6. LDA classification
#   7. Confusion matrices, ROC/AUC, LD scores, and variable importance
#   8. Robustness assessment by repeated train/test splits, permutation
#      testing, block-scaling comparison, and variable-importance stability
#
# Notes for reproducibility:
#   - The object Data_Fusion_3 may be loaded in the R environment, or the
#     dataset may be provided as a CSV file in data/processed/data_fusion_3.csv.
#   - The processed dataset must contain the response column "Grupos" and
#     all predictor variables listed below.
#   - Package versions should be reported using sessionInfo.txt or renv.lock.
# ============================================================


# ============================================================
# 1) SETTINGS AND PACKAGES
# ============================================================

required_packages <- c(
  "tidyverse",
  "caret",
  "MASS",
  "pROC",
  "rrcov",
  "ggrepel"
)

missing_packages <- setdiff(
  required_packages,
  rownames(installed.packages())
)

if (length(missing_packages) > 0) {
  stop(
    "The following packages are required but not installed: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

set.seed(123)

train_proportion <- 0.60
robpca_components <- 3
n_repeated_splits <- 100
n_permutations <- 500
top_n_variables <- 10

output_dir <- "outputs"
table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)


# ============================================================
# 2) DATA IMPORT
# ============================================================

input_file <- file.path("data", "processed", "data_fusion_3.csv")

if (file.exists(input_file)) {
  data_raw <- readr::read_csv(input_file, show_col_types = FALSE)
} else if (exists("Data_Fusion_3")) {
  data_raw <- Data_Fusion_3
} else {
  stop(
    "No input dataset was found. Provide a CSV file at ",
    input_file,
    " or load an object named Data_Fusion_3 in the R environment."
  )
}


# ============================================================
# 3) VARIABLE DEFINITIONS
# ============================================================

group_var <- "Grupos"

inorganic_vars <- c(
  "Zn", "Sr", "Ni", "Mn", "Fe", "Cu", "Co", "Cd", "Ba", "B"
)

organic_isotopic_vars <- c(
  "n-C18", "n-C19", "n-C21", "n-C24", "n-C29", "n-C30", "n-C31", "n-C33",
  "C29/C31", "ACL", "alcanos totais", "C18/total", "Delta C13", "%N", "%C"
)

fusion_vars <- c(inorganic_vars, organic_isotopic_vars)

variable_labels <- tibble::tibble(
  variable_original = fusion_vars,
  variable_label = c(
    "Zn", "Sr", "Ni", "Mn", "Fe", "Cu", "Co", "Cd", "Ba", "B",
    "n-C18", "n-C19", "n-C21", "n-C24", "n-C29", "n-C30", "n-C31", "n-C33",
    "C29/C31", "ACL", "Total alkanes", "C18/total", "delta13C",
    "%N", "%C"
  )
)


# ============================================================
# 4) HELPER FUNCTIONS
# ============================================================

check_required_columns <- function(data, required_columns) {
  missing_columns <- setdiff(required_columns, names(data))

  if (length(missing_columns) > 0) {
    stop(
      "The following required columns are missing from the dataset: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  invisible(TRUE)
}


to_numeric_dataframe <- function(data) {
  data <- as.data.frame(data)

  converted <- lapply(data, function(x) {
    if (is.factor(x)) {
      x <- as.character(x)
    }

    suppressWarnings(as.numeric(x))
  })

  converted <- as.data.frame(converted)
  names(converted) <- names(data)

  if (anyNA(converted) && !anyNA(data)) {
    warning(
      "Some values became NA after numeric conversion. ",
      "Check whether decimal separators or non-numeric symbols are present."
    )
  }

  converted
}


detect_zero_variance <- function(data) {
  data_numeric <- to_numeric_dataframe(data)

  sd_values <- apply(data_numeric, 2, sd, na.rm = TRUE)

  names(sd_values[is.na(sd_values) | sd_values == 0])
}


autoscale_data <- function(X) {
  X <- to_numeric_dataframe(X)

  means <- apply(X, 2, mean, na.rm = TRUE)
  sds <- apply(X, 2, sd, na.rm = TRUE)

  sds[is.na(sds) | sds == 0] <- 1

  X_scaled <- sweep(X, 2, means, "-")
  X_scaled <- sweep(X_scaled, 2, sds, "/")

  list(
    X_scaled = as.data.frame(X_scaled),
    means = means,
    sds = sds
  )
}


autoscale_train_test <- function(X_train, X_test) {
  X_train <- to_numeric_dataframe(X_train)
  X_test <- to_numeric_dataframe(X_test)

  zero_variance_vars <- detect_zero_variance(X_train)

  if (length(zero_variance_vars) > 0) {
    X_train <- X_train %>%
      dplyr::select(-dplyr::all_of(zero_variance_vars))

    X_test <- X_test %>%
      dplyr::select(-dplyr::all_of(zero_variance_vars))
  }

  means <- apply(X_train, 2, mean, na.rm = TRUE)
  sds <- apply(X_train, 2, sd, na.rm = TRUE)

  sds[is.na(sds) | sds == 0] <- 1

  X_train_scaled <- sweep(X_train, 2, means, "-")
  X_train_scaled <- sweep(X_train_scaled, 2, sds, "/")

  X_test_scaled <- sweep(X_test, 2, means, "-")
  X_test_scaled <- sweep(X_test_scaled, 2, sds, "/")

  list(
    train = as.data.frame(X_train_scaled),
    test = as.data.frame(X_test_scaled),
    means = means,
    sds = sds,
    removed_zero_variance_vars = zero_variance_vars
  )
}


make_model_names <- function(X_train, X_test) {
  original_names <- colnames(X_train)
  model_names <- make.names(original_names, unique = TRUE)

  colnames(X_train) <- model_names
  colnames(X_test) <- model_names

  variable_map <- tibble::tibble(
    variable_model = model_names,
    variable_original = original_names
  )

  list(
    X_train = X_train,
    X_test = X_test,
    variable_map = variable_map
  )
}


build_model_matrices <- function(train_data,
                                 test_data,
                                 model_type = c("inorganic", "organic_isotopic", "fusion"),
                                 block_scaling = TRUE,
                                 inorganic_vars_current,
                                 organic_isotopic_vars_current) {
  model_type <- match.arg(model_type)

  if (model_type == "inorganic") {
    X_train <- train_data %>%
      dplyr::select(dplyr::all_of(inorganic_vars_current))

    X_test <- test_data %>%
      dplyr::select(dplyr::all_of(inorganic_vars_current))

    scaled <- autoscale_train_test(X_train, X_test)

    return(list(
      X_train = scaled$train,
      X_test = scaled$test,
      removed_zero_variance_vars = scaled$removed_zero_variance_vars
    ))
  }

  if (model_type == "organic_isotopic") {
    X_train <- train_data %>%
      dplyr::select(dplyr::all_of(organic_isotopic_vars_current))

    X_test <- test_data %>%
      dplyr::select(dplyr::all_of(organic_isotopic_vars_current))

    scaled <- autoscale_train_test(X_train, X_test)

    return(list(
      X_train = scaled$train,
      X_test = scaled$test,
      removed_zero_variance_vars = scaled$removed_zero_variance_vars
    ))
  }

  X_inorg_train <- train_data %>%
    dplyr::select(dplyr::all_of(inorganic_vars_current))

  X_inorg_test <- test_data %>%
    dplyr::select(dplyr::all_of(inorganic_vars_current))

  X_org_train <- train_data %>%
    dplyr::select(dplyr::all_of(organic_isotopic_vars_current))

  X_org_test <- test_data %>%
    dplyr::select(dplyr::all_of(organic_isotopic_vars_current))

  scaled_inorg <- autoscale_train_test(X_inorg_train, X_inorg_test)
  scaled_org <- autoscale_train_test(X_org_train, X_org_test)

  X_inorg_train_scaled <- scaled_inorg$train
  X_inorg_test_scaled <- scaled_inorg$test

  X_org_train_scaled <- scaled_org$train
  X_org_test_scaled <- scaled_org$test

  if (block_scaling) {
    X_inorg_train_scaled <- X_inorg_train_scaled / sqrt(ncol(X_inorg_train_scaled))
    X_inorg_test_scaled <- X_inorg_test_scaled / sqrt(ncol(X_inorg_test_scaled))

    X_org_train_scaled <- X_org_train_scaled / sqrt(ncol(X_org_train_scaled))
    X_org_test_scaled <- X_org_test_scaled / sqrt(ncol(X_org_test_scaled))
  }

  list(
    X_train = as.data.frame(cbind(X_inorg_train_scaled, X_org_train_scaled)),
    X_test = as.data.frame(cbind(X_inorg_test_scaled, X_org_test_scaled)),
    removed_zero_variance_vars = c(
      scaled_inorg$removed_zero_variance_vars,
      scaled_org$removed_zero_variance_vars
    )
  )
}


calculate_class_metrics <- function(confusion_table) {
  classes <- rownames(confusion_table)

  metrics <- tibble::tibble(
    Class = classes,
    Sensitivity = NA_real_,
    Specificity = NA_real_,
    Precision = NA_real_,
    F1 = NA_real_
  )

  for (i in seq_along(classes)) {
    class_i <- classes[i]

    TP <- confusion_table[class_i, class_i]
    FN <- sum(confusion_table[class_i, ]) - TP
    FP <- sum(confusion_table[, class_i]) - TP
    TN <- sum(confusion_table) - TP - FN - FP

    sensitivity <- ifelse((TP + FN) == 0, NA_real_, TP / (TP + FN))
    specificity <- ifelse((TN + FP) == 0, NA_real_, TN / (TN + FP))
    precision <- ifelse((TP + FP) == 0, NA_real_, TP / (TP + FP))

    f1 <- ifelse(
      is.na(precision) | is.na(sensitivity) | (precision + sensitivity) == 0,
      NA_real_,
      2 * precision * sensitivity / (precision + sensitivity)
    )

    metrics$Sensitivity[i] <- sensitivity
    metrics$Specificity[i] <- specificity
    metrics$Precision[i] <- precision
    metrics$F1[i] <- f1
  }

  metrics
}


plot_confusion_matrix <- function(confusion_table, title_text) {
  cm_df <- as.data.frame(confusion_table)

  ggplot(cm_df, aes(x = Prediction, y = Reference, fill = Freq)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Freq), size = 5) +
    theme_bw(base_size = 14) +
    labs(
      title = title_text,
      x = "Predicted class",
      y = "Reference class",
      fill = "N"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}


summarize_lda_importance <- function(lda_model,
                                     variable_map,
                                     inorganic_vars_current,
                                     organic_isotopic_vars_current) {
  coefficients <- as.data.frame(lda_model$scaling)
  coefficients$variable_model <- rownames(coefficients)

  coefficients <- coefficients %>%
    dplyr::left_join(variable_map, by = "variable_model") %>%
    dplyr::mutate(
      block = dplyr::case_when(
        variable_original %in% inorganic_vars_current ~ "Inorganic",
        variable_original %in% organic_isotopic_vars_current ~ "Organic/isotopic",
        TRUE ~ "Unknown"
      )
    )

  ld_columns <- grep("^LD", names(coefficients), value = TRUE)

  coefficients <- coefficients %>%
    dplyr::mutate(
      global_importance = rowSums(abs(dplyr::across(dplyr::all_of(ld_columns))))
    ) %>%
    dplyr::arrange(dplyr::desc(global_importance))

  coefficients
}


fit_lda_once <- function(train_data,
                         test_data,
                         model_type = c("inorganic", "organic_isotopic", "fusion"),
                         block_scaling = TRUE,
                         seed = 123,
                         inorganic_vars_current,
                         organic_isotopic_vars_current) {
  model_type <- match.arg(model_type)

  set.seed(seed)

  y_train <- as.factor(train_data[[group_var]])
  y_test <- factor(test_data[[group_var]], levels = levels(y_train))

  matrices <- build_model_matrices(
    train_data = train_data,
    test_data = test_data,
    model_type = model_type,
    block_scaling = block_scaling,
    inorganic_vars_current = inorganic_vars_current,
    organic_isotopic_vars_current = organic_isotopic_vars_current
  )

  named <- make_model_names(matrices$X_train, matrices$X_test)

  X_train <- named$X_train
  X_test <- named$X_test
  variable_map <- named$variable_map

  lda_model <- MASS::lda(
    x = X_train,
    grouping = y_train
  )

  pred_train <- predict(lda_model, newdata = X_train)
  pred_test <- predict(lda_model, newdata = X_test)

  class_train <- factor(pred_train$class, levels = levels(y_train))
  class_test <- factor(pred_test$class, levels = levels(y_train))

  train_accuracy <- mean(class_train == y_train)
  test_accuracy <- mean(class_test == y_test)

  cm_train <- caret::confusionMatrix(
    data = class_train,
    reference = y_train
  )

  cm_test <- caret::confusionMatrix(
    data = class_test,
    reference = y_test
  )

  cm_train_table <- table(
    Reference = y_train,
    Prediction = class_train
  )

  cm_test_table <- table(
    Reference = y_test,
    Prediction = class_test
  )

  scores_train <- as.data.frame(pred_train$x)
  scores_test <- as.data.frame(pred_test$x)

  scores_train$Class <- y_train
  scores_train$Set <- "Training"

  scores_test$Class <- y_test
  scores_test$Set <- "Test"

  scores_all <- dplyr::bind_rows(scores_train, scores_test)

  posterior_train <- as.data.frame(pred_train$posterior)
  posterior_test <- as.data.frame(pred_test$posterior)

  eigenvalues <- lda_model$svd^2
  variance_explained <- eigenvalues / sum(eigenvalues) * 100

  variance_df <- tibble::tibble(
    LD = paste0("LD", seq_along(variance_explained)),
    Eigenvalue = eigenvalues,
    Explained_variance_percent = variance_explained,
    Cumulative_variance_percent = cumsum(variance_explained)
  )

  variable_importance <- summarize_lda_importance(
    lda_model = lda_model,
    variable_map = variable_map,
    inorganic_vars_current = inorganic_vars_current,
    organic_isotopic_vars_current = organic_isotopic_vars_current
  )

  list(
    model = lda_model,
    X_train = X_train,
    X_test = X_test,
    y_train = y_train,
    y_test = y_test,
    pred_train = pred_train,
    pred_test = pred_test,
    posterior_train = posterior_train,
    posterior_test = posterior_test,
    class_train = class_train,
    class_test = class_test,
    cm_train = cm_train,
    cm_test = cm_test,
    cm_train_table = cm_train_table,
    cm_test_table = cm_test_table,
    scores_all = scores_all,
    train_accuracy = train_accuracy,
    test_accuracy = test_accuracy,
    variance_df = variance_df,
    variable_importance = variable_importance,
    variable_map = variable_map,
    removed_zero_variance_vars = matrices$removed_zero_variance_vars
  )
}


run_single_split <- function(base_data,
                             model_type = c("inorganic", "organic_isotopic", "fusion"),
                             train_proportion = 0.60,
                             seed = 123,
                             block_scaling = TRUE,
                             permute_labels = FALSE,
                             inorganic_vars_current,
                             organic_isotopic_vars_current) {
  model_type <- match.arg(model_type)

  set.seed(seed)

  current_vars <- c(inorganic_vars_current, organic_isotopic_vars_current)

  base <- base_data %>%
    dplyr::select(dplyr::all_of(c(group_var, current_vars))) %>%
    tidyr::drop_na()

  base[[group_var]] <- as.factor(base[[group_var]])

  if (permute_labels) {
    base[[group_var]] <- sample(base[[group_var]])
    base[[group_var]] <- as.factor(base[[group_var]])
  }

  train_index <- caret::createDataPartition(
    base[[group_var]],
    p = train_proportion,
    list = FALSE
  )

  train_data <- base[train_index, , drop = FALSE]
  test_data <- base[-train_index, , drop = FALSE]

  result <- tryCatch(
    {
      fit <- fit_lda_once(
        train_data = train_data,
        test_data = test_data,
        model_type = model_type,
        block_scaling = block_scaling,
        seed = seed,
        inorganic_vars_current = inorganic_vars_current,
        organic_isotopic_vars_current = organic_isotopic_vars_current
      )

      variable_importance <- fit$variable_importance %>%
        dplyr::mutate(
          model_type = model_type,
          seed = seed,
          rank = dplyr::row_number()
        )

      list(
        success = TRUE,
        metrics = tibble::tibble(
          model_type = model_type,
          seed = seed,
          accuracy = fit$test_accuracy,
          n_train = nrow(train_data),
          n_test = nrow(test_data),
          block_scaling = block_scaling,
          permuted = permute_labels
        ),
        importance = variable_importance,
        predictions = tibble::tibble(
          observed = fit$y_test,
          predicted = fit$class_test,
          model_type = model_type,
          seed = seed
        )
      )
    },
    error = function(e) {
      message(
        "Error in split ", seed,
        " for model ", model_type,
        ": ", e$message
      )

      list(
        success = FALSE,
        metrics = tibble::tibble(
          model_type = model_type,
          seed = seed,
          accuracy = NA_real_,
          n_train = NA_integer_,
          n_test = NA_integer_,
          block_scaling = block_scaling,
          permuted = permute_labels
        ),
        importance = NULL,
        predictions = NULL
      )
    }
  )

  result
}


# ============================================================
# 5) DATA CHECKING AND PREPARATION
# ============================================================

check_required_columns(data_raw, c(group_var, fusion_vars))

data_model <- data_raw %>%
  dplyr::select(dplyr::all_of(c(group_var, fusion_vars))) %>%
  tidyr::drop_na()

data_model[[group_var]] <- as.factor(data_model[[group_var]])

cat("\nFull dataset dimensions before zero-variance filtering:\n")
print(dim(data_model))

cat("\nClass distribution before zero-variance filtering:\n")
print(table(data_model[[group_var]]))

zero_variance_vars_full <- detect_zero_variance(
  data_model %>%
    dplyr::select(dplyr::all_of(fusion_vars))
)

if (length(zero_variance_vars_full) > 0) {
  cat("\nVariables removed due to zero variance in the full dataset:\n")
  print(zero_variance_vars_full)

  data_model <- data_model %>%
    dplyr::select(-dplyr::all_of(zero_variance_vars_full))
}

inorganic_vars_current <- setdiff(inorganic_vars, zero_variance_vars_full)
organic_isotopic_vars_current <- setdiff(organic_isotopic_vars, zero_variance_vars_full)
fusion_vars_current <- c(inorganic_vars_current, organic_isotopic_vars_current)

readr::write_csv(
  tibble::tibble(removed_zero_variance_variable = zero_variance_vars_full),
  file.path(table_dir, "zero_variance_variables_removed.csv")
)


# ============================================================
# 6) ROBPCA-BASED OUTLIER SCREENING ON THE FUSED MATRIX
# ============================================================
# The ROBPCA model is fitted to the complete fused matrix after autoscaling
# within each analytical block and block scaling. This is used as an
# unsupervised multivariate screening step before the train/test split.

X_inorg_full <- data_model %>%
  dplyr::select(dplyr::all_of(inorganic_vars_current))

X_org_full <- data_model %>%
  dplyr::select(dplyr::all_of(organic_isotopic_vars_current))

scaled_inorg_full <- autoscale_data(X_inorg_full)
scaled_org_full <- autoscale_data(X_org_full)

X_inorg_full_scaled <- scaled_inorg_full$X_scaled / sqrt(ncol(scaled_inorg_full$X_scaled))
X_org_full_scaled <- scaled_org_full$X_scaled / sqrt(ncol(scaled_org_full$X_scaled))

X_fusion_full_scaled <- as.data.frame(
  cbind(X_inorg_full_scaled, X_org_full_scaled)
)

cat("\nDimensions of the fused matrix used for ROBPCA:\n")
print(dim(X_fusion_full_scaled))

robpca_fit <- rrcov::PcaHubert(
  X_fusion_full_scaled,
  k = robpca_components,
  scale = FALSE
)

robpca_flag <- robpca_fit@flag
outlier_flag <- robpca_flag == 0

cat(
  "\nNumber of multivariate outliers detected by ROBPCA:",
  sum(outlier_flag),
  "\n"
)

cat("\nIndices of ROBPCA outliers in the complete dataset:\n")
print(which(outlier_flag))

cat("\nDistribution of ROBPCA outliers by class:\n")
print(table(
  Class = data_model[[group_var]],
  Outlier = outlier_flag
))

outlier_table <- tibble::tibble(
  original_sample_index = seq_len(nrow(data_model)),
  class = data_model[[group_var]],
  robpca_flag = robpca_flag,
  outlier = outlier_flag
)

readr::write_csv(
  outlier_table,
  file.path(table_dir, "robpca_outliers_complete_dataset.csv")
)

png(
  filename = file.path(figure_dir, "robpca_diagnostic_plot.png"),
  width = 2400,
  height = 1800,
  res = 300
)
plot(robpca_fit)
dev.off()

data_clean <- data_model[!outlier_flag, , drop = FALSE]

cat("\nNumber of samples before ROBPCA outlier removal:", nrow(data_model), "\n")
cat("Number of samples after ROBPCA outlier removal:", nrow(data_clean), "\n")

cat("\nClass distribution after ROBPCA outlier removal:\n")
print(table(data_clean[[group_var]]))


# ============================================================
# 7) STRATIFIED TRAIN/TEST SPLIT
# ============================================================

set.seed(123)

train_index <- caret::createDataPartition(
  data_clean[[group_var]],
  p = train_proportion,
  list = FALSE
)

train_data <- data_clean[train_index, , drop = FALSE]
test_data <- data_clean[-train_index, , drop = FALSE]

cat("\nClass distribution in the training set:\n")
print(table(train_data[[group_var]]))

cat("\nClass distribution in the test set:\n")
print(table(test_data[[group_var]]))

split_table <- tibble::tibble(
  cleaned_sample_index = seq_len(nrow(data_clean)),
  class = data_clean[[group_var]],
  set = ifelse(seq_len(nrow(data_clean)) %in% train_index, "Training", "Test")
)

readr::write_csv(
  split_table,
  file.path(table_dir, "train_test_split_after_robpca.csv")
)


# ============================================================
# 8) MAIN LOW-LEVEL DATA FUSION LDA MODEL
# ============================================================

fusion_fit <- fit_lda_once(
  train_data = train_data,
  test_data = test_data,
  model_type = "fusion",
  block_scaling = TRUE,
  seed = 123,
  inorganic_vars_current = inorganic_vars_current,
  organic_isotopic_vars_current = organic_isotopic_vars_current
)

cat("\n==============================\n")
cat("LDA MODEL - LOW-LEVEL DATA FUSION\n")
cat("==============================\n")
print(fusion_fit$model)

cat("\nTraining confusion matrix:\n")
print(fusion_fit$cm_train)

cat("\nTest confusion matrix:\n")
print(fusion_fit$cm_test)

cat("\nTraining accuracy:", round(fusion_fit$train_accuracy * 100, 2), "%\n")
cat("Test accuracy:", round(fusion_fit$test_accuracy * 100, 2), "%\n")

prior_df <- tibble::tibble(
  Class = names(fusion_fit$model$prior),
  Prior_probability = as.numeric(fusion_fit$model$prior)
)

group_means_df <- as.data.frame(fusion_fit$model$means)
group_means_df$Class <- rownames(group_means_df)
group_means_df <- dplyr::relocate(group_means_df, Class)

class_metrics_test <- calculate_class_metrics(fusion_fit$cm_test_table)

readr::write_csv(
  as.data.frame(fusion_fit$cm_train_table),
  file.path(table_dir, "confusion_matrix_training_fusion.csv")
)

readr::write_csv(
  as.data.frame(fusion_fit$cm_test_table),
  file.path(table_dir, "confusion_matrix_test_fusion.csv")
)

readr::write_csv(
  fusion_fit$variable_importance,
  file.path(table_dir, "lda_coefficients_variable_importance_fusion.csv")
)

readr::write_csv(
  fusion_fit$variance_df,
  file.path(table_dir, "lda_explained_variance_fusion.csv")
)

readr::write_csv(
  class_metrics_test,
  file.path(table_dir, "test_class_metrics_fusion.csv")
)

readr::write_csv(
  fusion_fit$scores_all,
  file.path(table_dir, "lda_scores_training_test_fusion.csv")
)

readr::write_csv(
  prior_df,
  file.path(table_dir, "lda_prior_probabilities_fusion.csv")
)

readr::write_csv(
  group_means_df,
  file.path(table_dir, "lda_group_means_fusion.csv")
)


# ============================================================
# 9) MAIN MODEL FIGURES
# ============================================================

confusion_plot_test <- plot_confusion_matrix(
  fusion_fit$cm_test_table,
  "Confusion matrix - low-level data fusion LDA"
)

print(confusion_plot_test)

ggsave(
  filename = file.path(figure_dir, "confusion_matrix_test_fusion.png"),
  plot = confusion_plot_test,
  width = 7,
  height = 5.5,
  dpi = 300
)

if (all(c("LD1", "LD2") %in% colnames(fusion_fit$scores_all))) {
  var_exp <- fusion_fit$variance_df$Explained_variance_percent

  ld_plot <- ggplot(
    fusion_fit$scores_all,
    aes(x = LD1, y = LD2, color = Class, shape = Set)
  ) +
    geom_point(size = 3, alpha = 0.85) +
    theme_bw(base_size = 14) +
    labs(
      title = "LDA scores - low-level fused matrix",
      x = paste0("LD1 (", round(var_exp[1], 2), "%)"),
      y = paste0("LD2 (", round(var_exp[2], 2), "%)"),
      color = "Class",
      shape = "Set"
    )

  print(ld_plot)

  ggsave(
    filename = file.path(figure_dir, "lda_scores_ld1_ld2_fusion.png"),
    plot = ld_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
}

top_variables <- fusion_fit$variable_importance %>%
  dplyr::slice_head(n = min(20, nrow(.)))

variable_importance_plot <- ggplot(
  top_variables,
  aes(
    x = reorder(variable_original, global_importance),
    y = global_importance,
    fill = block
  )
) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 14) +
  labs(
    title = "Most important variables - low-level data fusion LDA",
    x = "Variable",
    y = "Global importance",
    fill = "Block"
  )

print(variable_importance_plot)

ggsave(
  filename = file.path(figure_dir, "variable_importance_top20_fusion.png"),
  plot = variable_importance_plot,
  width = 8,
  height = 6,
  dpi = 300
)


# ============================================================
# 10) ONE-VS-ALL ROC/AUC ANALYSIS FOR THE TEST SET
# ============================================================

classes <- levels(fusion_fit$y_train)

auc_values <- rep(NA_real_, length(classes))
names(auc_values) <- classes

roc_list <- list()

png(
  filename = file.path(figure_dir, "roc_one_vs_all_fusion.png"),
  width = 2400,
  height = 1800,
  res = 300
)

plot(
  NA,
  xlim = c(1, 0),
  ylim = c(0, 1),
  xlab = "Specificity",
  ylab = "Sensitivity",
  main = "One-vs-all ROC curves - low-level data fusion LDA"
)

abline(a = 1, b = -1, lty = 2, col = "gray")

for (i in seq_along(classes)) {
  class_i <- classes[i]

  if (!class_i %in% colnames(fusion_fit$posterior_test)) {
    warning("Class ", class_i, " is not present in posterior probabilities.")
    next
  }

  response_ova <- ifelse(fusion_fit$y_test == class_i, 1, 0)
  score_ova <- fusion_fit$posterior_test[[class_i]]

  if (length(unique(response_ova)) < 2) {
    warning(
      "Class ",
      class_i,
      " does not have both positive and negative samples in the test set."
    )
    next
  }

  roc_i <- pROC::roc(
    response = response_ova,
    predictor = score_ova,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )

  roc_list[[class_i]] <- roc_i
  auc_values[i] <- as.numeric(pROC::auc(roc_i))

  plot(
    roc_i,
    add = TRUE,
    legacy.axes = FALSE,
    print.auc = FALSE,
    col = i,
    lwd = 2
  )
}

legend(
  "bottomleft",
  legend = paste(names(auc_values), "- AUC =", round(auc_values, 3)),
  col = seq_along(classes),
  lwd = 2,
  cex = 0.8,
  bty = "n"
)

dev.off()

auc_df <- tibble::tibble(
  Class = names(auc_values),
  AUC_one_vs_all = as.numeric(auc_values)
)

cat("\nOne-vs-all AUC values for the test set:\n")
print(auc_df)

readr::write_csv(
  auc_df,
  file.path(table_dir, "auc_one_vs_all_fusion.csv")
)


# ============================================================
# 11) ROBUSTNESS ASSESSMENT: REPEATED TRAIN/TEST SPLITS
# ============================================================

models_to_compare <- c("inorganic", "organic_isotopic", "fusion")

repeated_results <- list()
counter <- 1

for (i in seq_len(n_repeated_splits)) {
  for (model_i in models_to_compare) {
    repeated_results[[counter]] <- run_single_split(
      base_data = data_clean,
      model_type = model_i,
      train_proportion = train_proportion,
      seed = 1000 + i,
      block_scaling = TRUE,
      permute_labels = FALSE,
      inorganic_vars_current = inorganic_vars_current,
      organic_isotopic_vars_current = organic_isotopic_vars_current
    )

    counter <- counter + 1
  }
}

repeated_metrics <- purrr::map_dfr(repeated_results, "metrics")
repeated_importance <- purrr::map_dfr(repeated_results, "importance")

robustness_summary <- repeated_metrics %>%
  dplyr::filter(!is.na(accuracy)) %>%
  dplyr::group_by(model_type) %>%
  dplyr::summarise(
    valid_repetitions = dplyr::n(),
    mean_accuracy = mean(accuracy),
    sd_accuracy = sd(accuracy),
    minimum_accuracy = min(accuracy),
    maximum_accuracy = max(accuracy),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    mean_accuracy_percent = 100 * mean_accuracy,
    sd_accuracy_percent = 100 * sd_accuracy,
    minimum_accuracy_percent = 100 * minimum_accuracy,
    maximum_accuracy_percent = 100 * maximum_accuracy
  )

cat("\nRepeated train/test validation summary:\n")
print(robustness_summary)

readr::write_csv(
  repeated_metrics,
  file.path(table_dir, "repeated_train_test_all_results.csv")
)

readr::write_csv(
  robustness_summary,
  file.path(table_dir, "repeated_train_test_summary.csv")
)

robustness_plot <- ggplot(
  repeated_metrics %>% dplyr::filter(!is.na(accuracy)),
  aes(x = model_type, y = accuracy * 100)
) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.5) +
  theme_bw(base_size = 14) +
  labs(
    title = "Repeated train/test validation",
    x = "Model",
    y = "Accuracy (%)"
  )

print(robustness_plot)

ggsave(
  filename = file.path(figure_dir, "repeated_train_test_accuracy_boxplot.png"),
  plot = robustness_plot,
  width = 8,
  height = 5,
  dpi = 300
)


# ============================================================
# 12) ROBUSTNESS ASSESSMENT: PERMUTATION TEST
# ============================================================

permutation_results <- list()

for (i in seq_len(n_permutations)) {
  permutation_results[[i]] <- run_single_split(
    base_data = data_clean,
    model_type = "fusion",
    train_proportion = train_proportion,
    seed = 5000 + i,
    block_scaling = TRUE,
    permute_labels = TRUE,
    inorganic_vars_current = inorganic_vars_current,
    organic_isotopic_vars_current = organic_isotopic_vars_current
  )
}

permutation_metrics <- purrr::map_dfr(permutation_results, "metrics") %>%
  dplyr::filter(!is.na(accuracy))

observed_fusion_accuracy <- robustness_summary %>%
  dplyr::filter(model_type == "fusion") %>%
  dplyr::pull(mean_accuracy)

empirical_p_value <- (
  sum(permutation_metrics$accuracy >= observed_fusion_accuracy) + 1
) / (nrow(permutation_metrics) + 1)

permutation_summary <- tibble::tibble(
  observed_fusion_mean_accuracy = observed_fusion_accuracy,
  mean_permuted_accuracy = mean(permutation_metrics$accuracy),
  sd_permuted_accuracy = sd(permutation_metrics$accuracy),
  minimum_permuted_accuracy = min(permutation_metrics$accuracy),
  maximum_permuted_accuracy = max(permutation_metrics$accuracy),
  empirical_p_value = empirical_p_value
)

cat("\nPermutation test summary:\n")
print(permutation_summary)

readr::write_csv(
  permutation_metrics,
  file.path(table_dir, "permutation_test_all_results.csv")
)

readr::write_csv(
  permutation_summary,
  file.path(table_dir, "permutation_test_summary.csv")
)

permutation_plot <- ggplot(permutation_metrics, aes(x = accuracy * 100)) +
  geom_histogram(bins = 30, color = "black") +
  geom_vline(
    xintercept = observed_fusion_accuracy * 100,
    linetype = "dashed",
    linewidth = 1
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = "Permutation test - low-level data fusion LDA",
    x = "Accuracy with permuted class labels (%)",
    y = "Frequency"
  )

print(permutation_plot)

ggsave(
  filename = file.path(figure_dir, "permutation_test_fusion.png"),
  plot = permutation_plot,
  width = 8,
  height = 5,
  dpi = 300
)


# ============================================================
# 13) ROBUSTNESS ASSESSMENT: EFFECT OF BLOCK SCALING
# ============================================================

block_scaling_results <- list()
counter <- 1

for (i in seq_len(n_repeated_splits)) {
  for (block_scaling_i in c(TRUE, FALSE)) {
    block_scaling_results[[counter]] <- run_single_split(
      base_data = data_clean,
      model_type = "fusion",
      train_proportion = train_proportion,
      seed = 8000 + i,
      block_scaling = block_scaling_i,
      permute_labels = FALSE,
      inorganic_vars_current = inorganic_vars_current,
      organic_isotopic_vars_current = organic_isotopic_vars_current
    )

    counter <- counter + 1
  }
}

block_scaling_metrics <- purrr::map_dfr(block_scaling_results, "metrics") %>%
  dplyr::filter(!is.na(accuracy)) %>%
  dplyr::mutate(
    strategy = ifelse(
      block_scaling,
      "With block scaling",
      "Without block scaling"
    )
  )

block_scaling_summary <- block_scaling_metrics %>%
  dplyr::group_by(strategy) %>%
  dplyr::summarise(
    valid_repetitions = dplyr::n(),
    mean_accuracy = mean(accuracy),
    sd_accuracy = sd(accuracy),
    minimum_accuracy = min(accuracy),
    maximum_accuracy = max(accuracy),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    mean_accuracy_percent = 100 * mean_accuracy,
    sd_accuracy_percent = 100 * sd_accuracy,
    minimum_accuracy_percent = 100 * minimum_accuracy,
    maximum_accuracy_percent = 100 * maximum_accuracy
  )

cat("\nBlock-scaling comparison:\n")
print(block_scaling_summary)

readr::write_csv(
  block_scaling_metrics,
  file.path(table_dir, "block_scaling_comparison_all_results.csv")
)

readr::write_csv(
  block_scaling_summary,
  file.path(table_dir, "block_scaling_comparison_summary.csv")
)

block_scaling_plot <- ggplot(
  block_scaling_metrics,
  aes(x = strategy, y = accuracy * 100)
) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.5) +
  theme_bw(base_size = 14) +
  labs(
    title = "Effect of block scaling on low-level data fusion",
    x = "",
    y = "Accuracy (%)"
  )

print(block_scaling_plot)

ggsave(
  filename = file.path(figure_dir, "block_scaling_comparison.png"),
  plot = block_scaling_plot,
  width = 8,
  height = 5,
  dpi = 300
)


# ============================================================
# 14) ROBUSTNESS ASSESSMENT: VARIABLE-IMPORTANCE STABILITY
# ============================================================

valid_fusion_seeds <- repeated_importance %>%
  dplyr::filter(model_type == "fusion") %>%
  dplyr::distinct(seed) %>%
  nrow()

variable_stability <- repeated_importance %>%
  dplyr::filter(model_type == "fusion") %>%
  dplyr::group_by(seed) %>%
  dplyr::arrange(rank, .by_group = TRUE) %>%
  dplyr::slice_head(n = top_n_variables) %>%
  dplyr::ungroup() %>%
  dplyr::count(variable_original, name = "frequency_in_top_n") %>%
  dplyr::mutate(
    frequency_percent = 100 * frequency_in_top_n / valid_fusion_seeds
  ) %>%
  dplyr::arrange(dplyr::desc(frequency_percent))

cat("\nVariable-importance stability for the fused model:\n")
print(variable_stability)

readr::write_csv(
  variable_stability,
  file.path(
    table_dir,
    paste0("variable_importance_stability_top", top_n_variables, ".csv")
  )
)

variable_stability_plot <- ggplot(
  variable_stability,
  aes(
    x = reorder(variable_original, frequency_percent),
    y = frequency_percent
  )
) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 14) +
  labs(
    title = paste0("Variable-importance stability - top ", top_n_variables),
    x = "Variable",
    y = paste0("Frequency in the top ", top_n_variables, " variables (%)")
  )

print(variable_stability_plot)

ggsave(
  filename = file.path(
    figure_dir,
    paste0("variable_importance_stability_top", top_n_variables, ".png")
  ),
  plot = variable_stability_plot,
  width = 8,
  height = 6,
  dpi = 300
)


# ============================================================
# 15) FINAL ARTICLE-READY TABLE
# ============================================================

article_robustness_table <- robustness_summary %>%
  dplyr::transmute(
    Model = dplyr::case_when(
      model_type == "inorganic" ~ "Inorganic LDA",
      model_type == "organic_isotopic" ~ "Organic/isotopic LDA",
      model_type == "fusion" ~ "Low-level data fusion LDA",
      TRUE ~ model_type
    ),
    `Mean accuracy (%)` = round(mean_accuracy_percent, 2),
    `SD (%)` = round(sd_accuracy_percent, 2),
    `Minimum accuracy (%)` = round(minimum_accuracy_percent, 2),
    `Maximum accuracy (%)` = round(maximum_accuracy_percent, 2)
  )

cat("\nArticle-ready robustness table:\n")
print(article_robustness_table)

readr::write_csv(
  article_robustness_table,
  file.path(table_dir, "article_ready_robustness_table.csv")
)


# ============================================================
# 16) COMPUTATIONAL ENVIRONMENT
# ============================================================

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)

cat("\nLow-level data fusion analysis completed successfully.\n")
cat("Tables saved in: ", table_dir, "\n", sep = "")
cat("Figures saved in: ", figure_dir, "\n", sep = "")
