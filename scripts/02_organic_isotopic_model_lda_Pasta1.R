# ============================================================
# ORGANIC/ISOTOPIC LDA MODEL FOR COMMERCIAL CIGARETTE CLASSIFICATION
# Supplementary reference material for manuscript submission
# ============================================================
#
# This script reproduces the organic/isotopic chemometric workflow used for
# commercial cigarette classification.
#
# Main steps:
#   1. Load packages and input data
#   2. Select organic/isotopic descriptors
#   3. Apply log1p transformation to selected n-alkane variables
#   4. Remove incomplete rows and zero-variance variables
#   5. Detect multivariate outliers by group using ROBPCA
#   6. Fit and validate LDA models
#   7. Export confusion matrices, LDA scores, ROC/AUC results,
#      variable-importance tables, robustness summaries, and figures
#
# Expected input:
#   Option 1: an object named Pasta1 already loaded in R
#   Option 2: a CSV file named pasta1.csv in data/processed/ or data/
#   Option 3: a CSV file named data_fusion_3.csv in data/processed/ or data/
#             if the same processed dataset is used for the article
#
# Expected group column:
#   - Grupos or grupos
#
# If no group column is present and the dataset contains 96 rows, the script
# can assign the group labels according to the original sample order used in
# the manuscript dataset.
#
# ============================================================


# ------------------------------------------------------------
# 1. Global settings
# ------------------------------------------------------------

set.seed(123)

output_dir <- "outputs"
tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

n_repeated_splits <- 100
n_permutations <- 500
training_fraction <- 0.60


# ------------------------------------------------------------
# 2. Packages
# ------------------------------------------------------------

required_packages <- c(
  "tidyverse",
  "caret",
  "MASS",
  "pROC",
  "rrcov",
  "ggplot2",
  "ggrepel"
)

packages_to_install <- required_packages[
  !(required_packages %in% rownames(installed.packages()))
]

if (length(packages_to_install) > 0) {
  install.packages(packages_to_install)
}

invisible(lapply(required_packages, library, character.only = TRUE))


# ------------------------------------------------------------
# 3. Input data
# ------------------------------------------------------------

if (exists("Pasta1")) {
  raw_data <- Pasta1
} else if (file.exists("data/processed/pasta1.csv")) {
  raw_data <- read.csv(
    "data/processed/pasta1.csv",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
} else if (file.exists("data/pasta1.csv")) {
  raw_data <- read.csv(
    "data/pasta1.csv",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
} else if (file.exists("data/processed/data_fusion_3.csv")) {
  raw_data <- read.csv(
    "data/processed/data_fusion_3.csv",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
} else if (file.exists("data/data_fusion_3.csv")) {
  raw_data <- read.csv(
    "data/data_fusion_3.csv",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
} else {
  stop(
    "Input data not found. Please load Pasta1 in R, ",
    "or place pasta1.csv in data/processed/ or data/."
  )
}

raw_data <- as.data.frame(raw_data, check.names = FALSE)


# ------------------------------------------------------------
# 4. Organic/isotopic descriptors
# ------------------------------------------------------------

organic_isotopic_vars <- c(
  "n-C18",
  "n-C19",
  "n-C21",
  "n-C24",
  "n-C29",
  "n-C30",
  "n-C31",
  "n-C33",
  "C29/C31",
  "ACL",
  "alcanos totais",
  "C18/total",
  "Delta C13",
  "%N",
  "%C"
)

log1p_vars <- c(
  "n-C18",
  "n-C19",
  "n-C21",
  "n-C24",
  "n-C29",
  "n-C30",
  "n-C31",
  "n-C33",
  "alcanos totais"
)

missing_vars <- setdiff(organic_isotopic_vars, names(raw_data))

if (length(missing_vars) > 0) {
  stop(
    "The following organic/isotopic variables are missing from the dataset: ",
    paste(missing_vars, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 5. Group handling
# ------------------------------------------------------------

multiclass_labels_original_order <- c(
  rep("Smuggled", 9),
  rep("Multinational", 21),
  rep("Smuggled", 6),
  rep("Regional", 3),
  rep("Smuggled", 3),
  rep("Multinational", 9),
  rep("Smuggled", 3),
  rep("Counterfeit", 3),
  rep("No_Health_Registration", 6),
  rep("Regional", 15),
  rep("No_Health_Registration", 3),
  rep("Regional", 6),
  rep("Counterfeit", 6),
  rep("Regional", 3)
)

standardize_group_labels <- function(x) {
  x <- as.character(x)

  x <- dplyr::case_when(
    x %in% c("Contrabandeado", "Smuggled") ~ "Smuggled",
    x %in% c("Falsificado", "Counterfeit") ~ "Counterfeit",
    x %in% c(
      "Sem_Registro_Sanitário",
      "Sem Registro Sanitário",
      "Sem_Registro_Sanitario",
      "Sem Registro Sanitario",
      "No_Health_Registration"
    ) ~ "No_Health_Registration",
    x %in% c("Regional", "Nacional") ~ "Regional",
    x %in% c("Multinacional", "Multinational") ~ "Multinational",
    x %in% c("Outros", "Others", "Other") ~ "Others",
    TRUE ~ x
  )

  factor(x)
}

if ("Grupos" %in% names(raw_data)) {
  raw_data$Group <- standardize_group_labels(raw_data$Grupos)
} else if ("grupos" %in% names(raw_data)) {
  raw_data$Group <- standardize_group_labels(raw_data$grupos)
} else if (nrow(raw_data) == length(multiclass_labels_original_order)) {
  raw_data$Group <- factor(multiclass_labels_original_order)
  warning(
    "No group column was found. Group labels were assigned according to ",
    "the original 96-sample order used in the manuscript dataset."
  )
} else {
  stop(
    "No group column was found. Please include a column named 'Grupos' or 'grupos'."
  )
}

cat("\nInitial class distribution:\n")
print(table(raw_data$Group))


# ------------------------------------------------------------
# 6. Build organic/isotopic modelling dataset
# ------------------------------------------------------------

model_data <- raw_data %>%
  dplyr::select(
    Group,
    dplyr::all_of(organic_isotopic_vars)
  )

for (v in organic_isotopic_vars) {
  model_data[[v]] <- as.numeric(model_data[[v]])
}

for (v in intersect(log1p_vars, names(model_data))) {
  model_data[[v]] <- log1p(model_data[[v]])
}

model_data <- model_data %>%
  tidyr::drop_na()

cat("\nDataset dimensions after log1p transformation and missing-value removal:\n")
print(dim(model_data))

cat("\nClass distribution after missing-value removal:\n")
print(table(model_data$Group))


# ------------------------------------------------------------
# 7. Remove zero-variance variables
# ------------------------------------------------------------

predictor_matrix <- model_data %>%
  dplyr::select(dplyr::all_of(organic_isotopic_vars)) %>%
  as.data.frame()

sd_vars <- apply(predictor_matrix, 2, sd, na.rm = TRUE)
zero_variance_vars <- names(sd_vars[sd_vars == 0 | is.na(sd_vars)])

if (length(zero_variance_vars) > 0) {
  cat("\nVariables removed due to zero variance:\n")
  print(zero_variance_vars)

  organic_isotopic_vars <- setdiff(organic_isotopic_vars, zero_variance_vars)

  model_data <- model_data %>%
    dplyr::select(Group, dplyr::all_of(organic_isotopic_vars))
}

cat("\nFinal organic/isotopic descriptors used in the models:\n")
print(organic_isotopic_vars)


# ------------------------------------------------------------
# 8. Descriptive summary
# ------------------------------------------------------------

descriptive_summary <- model_data %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(organic_isotopic_vars),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  dplyr::group_by(Group, Variable) %>%
  dplyr::summarise(
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    Median = median(Value, na.rm = TRUE),
    Minimum = min(Value, na.rm = TRUE),
    Maximum = max(Value, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  )

write.csv(
  descriptive_summary,
  file.path(tables_dir, "organic_isotopic_descriptive_summary.csv"),
  row.names = FALSE
)

descriptive_plot <- descriptive_summary %>%
  dplyr::mutate(
    SD = ifelse(is.na(SD), 0, SD)
  ) %>%
  ggplot2::ggplot(
    ggplot2::aes(x = Group, y = Mean, colour = Group)
  ) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = Mean - SD, ymax = Mean + SD),
    width = 0.15
  ) +
  ggplot2::facet_wrap(~ Variable, scales = "free_y", ncol = 3) +
  ggplot2::labs(
    x = NULL,
    y = "Mean +/- SD",
    colour = "Class"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    legend.position = "bottom",
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )

ggplot2::ggsave(
  file.path(figures_dir, "organic_isotopic_descriptive_summary.png"),
  descriptive_plot,
  width = 12,
  height = 8,
  dpi = 300
)


# ------------------------------------------------------------
# 9. ROBPCA-based outlier detection by group
# ------------------------------------------------------------

safe_slot <- function(object, slot_name, expected_length) {
  if (slot_name %in% slotNames(object)) {
    value <- slot(object, slot_name)
    value <- as.numeric(value)

    if (length(value) == expected_length) {
      return(value)
    }
  }

  rep(NA_real_, expected_length)
}

detect_robpca_outliers_single_group <- function(group_data, vars) {
  out <- group_data

  out$robpca_flag <- NA_integer_
  out$robpca_score_distance <- NA_real_
  out$robpca_orthogonal_distance <- NA_real_
  out$robpca_cutoff_score_distance <- NA_real_
  out$robpca_cutoff_orthogonal_distance <- NA_real_
  out$robpca_outlier <- NA
  out$robpca_note <- NA_character_

  X <- out[, vars, drop = FALSE]
  X <- as.data.frame(lapply(X, as.numeric))

  keep_vars <- sapply(X, function(z) {
    sum(!is.na(z)) > 1 && sd(z, na.rm = TRUE) > 0
  })

  X <- X[, keep_vars, drop = FALSE]

  if (ncol(X) < 2) {
    out$robpca_note <- "ROBPCA not applied: fewer than two variables with variability."
    return(out)
  }

  complete_rows <- complete.cases(X)

  if (sum(complete_rows) < 4) {
    out$robpca_note <- "ROBPCA not applied: fewer than four complete observations."
    return(out)
  }

  X_complete <- as.matrix(X[complete_rows, , drop = FALSE])
  X_scaled <- scale(X_complete)

  n_complete <- nrow(X_scaled)
  p_complete <- ncol(X_scaled)
  k_components <- min(2, p_complete - 1, n_complete - 2)

  if (k_components < 1) {
    out$robpca_note <- "ROBPCA not applied: insufficient dimensionality."
    return(out)
  }

  fit <- tryCatch(
    rrcov::PcaHubert(
      X_scaled,
      k = k_components,
      scale = FALSE
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    out$robpca_note <- paste("ROBPCA error:", fit$message)
    return(out)
  }

  flag <- as.integer(fit@flag)

  out$robpca_flag[complete_rows] <- flag
  out$robpca_score_distance[complete_rows] <- safe_slot(fit, "sd", n_complete)
  out$robpca_orthogonal_distance[complete_rows] <- safe_slot(fit, "od", n_complete)

  if ("cutoff.sd" %in% slotNames(fit)) {
    out$robpca_cutoff_score_distance[complete_rows] <- as.numeric(fit@cutoff.sd)[1]
  }

  if ("cutoff.od" %in% slotNames(fit)) {
    out$robpca_cutoff_orthogonal_distance[complete_rows] <- as.numeric(fit@cutoff.od)[1]
  }

  out$robpca_outlier[complete_rows] <- flag == 0
  out$robpca_note <- paste0(
    "ROBPCA applied using ",
    p_complete,
    " variables and ",
    k_components,
    " components."
  )

  out
}

robpca_results <- model_data %>%
  dplyr::group_by(Group) %>%
  dplyr::group_split(.keep = TRUE) %>%
  purrr::map_dfr(
    detect_robpca_outliers_single_group,
    vars = organic_isotopic_vars
  )

outlier_table <- robpca_results %>%
  dplyr::mutate(Sample_ID = dplyr::row_number()) %>%
  dplyr::select(
    Sample_ID,
    Group,
    robpca_flag,
    robpca_score_distance,
    robpca_orthogonal_distance,
    robpca_cutoff_score_distance,
    robpca_cutoff_orthogonal_distance,
    robpca_outlier,
    robpca_note
  )

write.csv(
  outlier_table,
  file.path(tables_dir, "organic_isotopic_robpca_outliers.csv"),
  row.names = FALSE
)

cat("\nROBPCA outliers by group:\n")
print(
  outlier_table %>%
    dplyr::filter(robpca_outlier == TRUE) %>%
    dplyr::count(Group, name = "Outliers")
)

clean_data <- robpca_results %>%
  dplyr::filter(is.na(robpca_outlier) | robpca_outlier == FALSE) %>%
  dplyr::select(Group, dplyr::all_of(organic_isotopic_vars))

clean_data$Group <- droplevels(as.factor(clean_data$Group))

cat("\nDataset dimensions before outlier removal:\n")
print(dim(model_data))

cat("\nDataset dimensions after outlier removal:\n")
print(dim(clean_data))

cat("\nClass distribution after outlier removal:\n")
print(table(clean_data$Group))

write.csv(
  clean_data,
  file.path(tables_dir, "organic_isotopic_clean_dataset_after_robpca.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. Utility functions for LDA modelling
# ------------------------------------------------------------

autoscale_train_test <- function(X_train, X_test) {
  X_train <- as.data.frame(lapply(X_train, as.numeric))
  X_test <- as.data.frame(lapply(X_test, as.numeric))

  train_means <- apply(X_train, 2, mean, na.rm = TRUE)
  train_sds <- apply(X_train, 2, sd, na.rm = TRUE)

  train_sds[train_sds == 0 | is.na(train_sds)] <- 1

  X_train_scaled <- sweep(X_train, 2, train_means, "-")
  X_train_scaled <- sweep(X_train_scaled, 2, train_sds, "/")

  X_test_scaled <- sweep(X_test, 2, train_means, "-")
  X_test_scaled <- sweep(X_test_scaled, 2, train_sds, "/")

  list(
    train = as.data.frame(X_train_scaled),
    test = as.data.frame(X_test_scaled),
    means = train_means,
    sds = train_sds
  )
}

calculate_class_metrics <- function(confusion_table) {
  classes <- rownames(confusion_table)

  metrics <- data.frame(
    Class = classes,
    Sensitivity = NA_real_,
    Specificity = NA_real_,
    Precision = NA_real_,
    F1 = NA_real_
  )

  for (i in seq_along(classes)) {
    current_class <- classes[i]

    TP <- confusion_table[current_class, current_class]
    FN <- sum(confusion_table[current_class, ]) - TP
    FP <- sum(confusion_table[, current_class]) - TP
    TN <- sum(confusion_table) - TP - FN - FP

    sensitivity <- ifelse((TP + FN) == 0, NA, TP / (TP + FN))
    specificity <- ifelse((TN + FP) == 0, NA, TN / (TN + FP))
    precision <- ifelse((TP + FP) == 0, NA, TP / (TP + FP))
    f1 <- ifelse(
      is.na(precision) | is.na(sensitivity) | (precision + sensitivity) == 0,
      NA,
      2 * precision * sensitivity / (precision + sensitivity)
    )

    metrics$Sensitivity[i] <- sensitivity
    metrics$Specificity[i] <- specificity
    metrics$Precision[i] <- precision
    metrics$F1[i] <- f1
  }

  metrics
}

calculate_roc_auc_one_vs_all <- function(observed, posterior_probabilities) {
  observed <- as.factor(observed)
  classes <- levels(observed)

  auc_values <- numeric(length(classes))
  names(auc_values) <- classes

  roc_list <- list()

  for (current_class in classes) {
    if (!current_class %in% colnames(posterior_probabilities)) {
      auc_values[current_class] <- NA_real_
      next
    }

    binary_response <- ifelse(observed == current_class, 1, 0)

    if (length(unique(binary_response)) < 2) {
      auc_values[current_class] <- NA_real_
      next
    }

    roc_object <- pROC::roc(
      response = binary_response,
      predictor = posterior_probabilities[, current_class],
      quiet = TRUE
    )

    roc_list[[current_class]] <- roc_object
    auc_values[current_class] <- as.numeric(pROC::auc(roc_object))
  }

  list(
    auc_table = data.frame(
      Class = names(auc_values),
      AUC_one_vs_all = as.numeric(auc_values),
      row.names = NULL
    ),
    roc_list = roc_list
  )
}

fit_validate_lda <- function(data,
                             outcome_col,
                             predictor_vars,
                             model_name,
                             output_prefix,
                             p_train = 0.60,
                             seed = 123) {
  set.seed(seed)

  modelling_data <- data %>%
    dplyr::select(
      dplyr::all_of(outcome_col),
      dplyr::all_of(predictor_vars)
    ) %>%
    tidyr::drop_na()

  colnames(modelling_data)[1] <- "Outcome"
  modelling_data$Outcome <- droplevels(as.factor(modelling_data$Outcome))

  if (nlevels(modelling_data$Outcome) < 2) {
    stop("The outcome variable must contain at least two classes.")
  }

  split_index <- caret::createDataPartition(
    modelling_data$Outcome,
    p = p_train,
    list = FALSE
  )

  training_data <- modelling_data[split_index, , drop = FALSE]
  test_data <- modelling_data[-split_index, , drop = FALSE]

  y_train <- droplevels(as.factor(training_data$Outcome))
  y_test <- factor(test_data$Outcome, levels = levels(y_train))

  X_train_raw <- training_data %>%
    dplyr::select(dplyr::all_of(predictor_vars))

  X_test_raw <- test_data %>%
    dplyr::select(dplyr::all_of(predictor_vars))

  scaled_data <- autoscale_train_test(X_train_raw, X_test_raw)

  X_train <- scaled_data$train
  X_test <- scaled_data$test

  original_names <- colnames(X_train)
  model_names <- make.names(original_names, unique = TRUE)

  colnames(X_train) <- model_names
  colnames(X_test) <- model_names

  variable_map <- data.frame(
    Model_variable = model_names,
    Original_variable = original_names
  )

  lda_model <- MASS::lda(
    x = X_train,
    grouping = y_train
  )

  train_prediction <- predict(lda_model, newdata = X_train)
  test_prediction <- predict(lda_model, newdata = X_test)

  train_classes <- factor(train_prediction$class, levels = levels(y_train))
  test_classes <- factor(test_prediction$class, levels = levels(y_train))

  train_probabilities <- as.data.frame(train_prediction$posterior)
  test_probabilities <- as.data.frame(test_prediction$posterior)

  train_scores <- as.data.frame(train_prediction$x)
  test_scores <- as.data.frame(test_prediction$x)

  train_confusion_table <- table(
    Observed = y_train,
    Predicted = train_classes
  )

  test_confusion_table <- table(
    Observed = y_test,
    Predicted = test_classes
  )

  train_accuracy <- mean(train_classes == y_train)
  test_accuracy <- mean(test_classes == y_test)

  class_metrics <- calculate_class_metrics(test_confusion_table)

  eigenvalues <- lda_model$svd^2
  explained_variance <- eigenvalues / sum(eigenvalues) * 100

  variance_table <- data.frame(
    LD = paste0("LD", seq_along(explained_variance)),
    Eigenvalue = eigenvalues,
    Explained_variance_percent = explained_variance,
    Cumulative_variance_percent = cumsum(explained_variance)
  )

  coefficients <- as.data.frame(lda_model$scaling)
  coefficients$Model_variable <- rownames(coefficients)

  coefficients <- coefficients %>%
    dplyr::left_join(variable_map, by = "Model_variable")

  ld_columns <- grep("^LD", names(coefficients), value = TRUE)

  coefficients$Global_importance <- rowSums(
    abs(coefficients[, ld_columns, drop = FALSE])
  )

  coefficients <- coefficients %>%
    dplyr::arrange(dplyr::desc(Global_importance))

  roc_results <- calculate_roc_auc_one_vs_all(
    observed = y_test,
    posterior_probabilities = test_probabilities
  )

  auc_table <- roc_results$auc_table

  write.csv(
    as.data.frame(train_confusion_table),
    file.path(tables_dir, paste0(output_prefix, "_training_confusion_matrix.csv")),
    row.names = FALSE
  )

  write.csv(
    as.data.frame(test_confusion_table),
    file.path(tables_dir, paste0(output_prefix, "_test_confusion_matrix.csv")),
    row.names = FALSE
  )

  write.csv(
    data.frame(
      Model = model_name,
      Training_accuracy = train_accuracy,
      Test_accuracy = test_accuracy,
      Training_n = nrow(training_data),
      Test_n = nrow(test_data)
    ),
    file.path(tables_dir, paste0(output_prefix, "_accuracy_summary.csv")),
    row.names = FALSE
  )

  write.csv(
    class_metrics,
    file.path(tables_dir, paste0(output_prefix, "_test_class_metrics.csv")),
    row.names = FALSE
  )

  write.csv(
    coefficients,
    file.path(tables_dir, paste0(output_prefix, "_variable_importance.csv")),
    row.names = FALSE
  )

  write.csv(
    variance_table,
    file.path(tables_dir, paste0(output_prefix, "_lda_explained_variance.csv")),
    row.names = FALSE
  )

  write.csv(
    auc_table,
    file.path(tables_dir, paste0(output_prefix, "_one_vs_all_auc.csv")),
    row.names = FALSE
  )

  test_confusion_df <- as.data.frame(test_confusion_table)

  confusion_plot <- ggplot2::ggplot(
    test_confusion_df,
    ggplot2::aes(x = Predicted, y = Observed, fill = Freq)
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = Freq), size = 5) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::labs(
      title = paste0(model_name, " - test confusion matrix"),
      x = "Predicted class",
      y = "Observed class",
      fill = "n"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  ggplot2::ggsave(
    file.path(figures_dir, paste0(output_prefix, "_test_confusion_matrix.png")),
    confusion_plot,
    width = 7,
    height = 6,
    dpi = 300
  )

  train_scores$Group <- y_train
  train_scores$Set <- "Training"

  test_scores$Group <- y_test
  test_scores$Set <- "Test"

  score_table <- dplyr::bind_rows(train_scores, test_scores)

  write.csv(
    score_table,
    file.path(tables_dir, paste0(output_prefix, "_lda_scores.csv")),
    row.names = FALSE
  )

  if ("LD2" %in% colnames(score_table)) {
    score_plot <- ggplot2::ggplot(
      score_table,
      ggplot2::aes(x = LD1, y = LD2, colour = Group, shape = Set)
    ) +
      ggplot2::geom_point(size = 3, alpha = 0.85) +
      ggplot2::theme_bw(base_size = 14) +
      ggplot2::labs(
        title = paste0(model_name, " - LDA scores"),
        x = paste0("LD1 (", round(explained_variance[1], 2), "%)"),
        y = paste0("LD2 (", round(explained_variance[2], 2), "%)"),
        colour = "Class",
        shape = "Set"
      )
  } else {
    score_plot <- ggplot2::ggplot(
      score_table,
      ggplot2::aes(x = LD1, y = 0, colour = Group, shape = Set)
    ) +
      ggplot2::geom_point(size = 3, alpha = 0.85) +
      ggplot2::theme_bw(base_size = 14) +
      ggplot2::labs(
        title = paste0(model_name, " - LDA scores"),
        x = "LD1",
        y = NULL,
        colour = "Class",
        shape = "Set"
      ) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
  }

  ggplot2::ggsave(
    file.path(figures_dir, paste0(output_prefix, "_lda_scores.png")),
    score_plot,
    width = 8,
    height = 6,
    dpi = 300
  )

  top_variables <- coefficients %>%
    dplyr::slice_head(n = min(20, dplyr::n()))

  importance_plot <- ggplot2::ggplot(
    top_variables,
    ggplot2::aes(
      x = reorder(Original_variable, Global_importance),
      y = Global_importance
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::labs(
      title = paste0(model_name, " - variable importance"),
      x = "Variable",
      y = "Global importance"
    )

  ggplot2::ggsave(
    file.path(figures_dir, paste0(output_prefix, "_variable_importance.png")),
    importance_plot,
    width = 8,
    height = 6,
    dpi = 300
  )

  roc_png_path <- file.path(figures_dir, paste0(output_prefix, "_roc_one_vs_all.png"))

  png(roc_png_path, width = 2400, height = 1800, res = 300)
  plot(
    NA,
    xlim = c(1, 0),
    ylim = c(0, 1),
    xlab = "Specificity",
    ylab = "Sensitivity",
    main = paste0(model_name, " - one-vs-all ROC curves")
  )
  abline(a = 1, b = -1, lty = 2, col = "gray")

  valid_roc_names <- names(roc_results$roc_list)

  if (length(valid_roc_names) > 0) {
    for (current_class in valid_roc_names) {
      plot(
        roc_results$roc_list[[current_class]],
        add = TRUE,
        print.auc = FALSE
      )
    }

    legend(
      "bottomleft",
      legend = paste0(
        auc_table$Class,
        " - AUC = ",
        round(auc_table$AUC_one_vs_all, 3)
      ),
      cex = 0.8,
      bty = "n"
    )
  }
  dev.off()

  list(
    model = lda_model,
    train_confusion_table = train_confusion_table,
    test_confusion_table = test_confusion_table,
    class_metrics = class_metrics,
    coefficients = coefficients,
    variance_table = variance_table,
    auc_table = auc_table,
    train_accuracy = train_accuracy,
    test_accuracy = test_accuracy,
    score_table = score_table
  )
}


# ------------------------------------------------------------
# 11. Fit multiclass organic/isotopic LDA model
# ------------------------------------------------------------

multiclass_results <- NULL

if (nlevels(clean_data$Group) > 2) {
  cat("\nFitting multiclass organic/isotopic LDA model...\n")

  multiclass_results <- fit_validate_lda(
    data = clean_data,
    outcome_col = "Group",
    predictor_vars = organic_isotopic_vars,
    model_name = "Organic/isotopic LDA model",
    output_prefix = "organic_isotopic_multiclass_lda",
    p_train = training_fraction,
    seed = 123
  )

  cat("\nMulticlass organic/isotopic LDA test accuracy:\n")
  print(multiclass_results$test_accuracy)
} else {
  cat("\nMulticlass LDA skipped because the dataset contains only two classes.\n")
}


# ------------------------------------------------------------
# 12. Fit binary organic/isotopic LDA model
# ------------------------------------------------------------

binary_data <- clean_data %>%
  dplyr::mutate(
    Binary_Group = dplyr::if_else(
      Group == "Multinational",
      "Multinational",
      "Others"
    )
  )

binary_data$Binary_Group <- factor(
  binary_data$Binary_Group,
  levels = c("Others", "Multinational")
)

binary_results <- NULL

if (nlevels(droplevels(binary_data$Binary_Group)) == 2) {
  cat("\nFitting binary organic/isotopic LDA model: Multinational vs Others...\n")

  binary_results <- fit_validate_lda(
    data = binary_data,
    outcome_col = "Binary_Group",
    predictor_vars = organic_isotopic_vars,
    model_name = "Organic/isotopic LDA model - Multinational vs Others",
    output_prefix = "organic_isotopic_binary_lda",
    p_train = training_fraction,
    seed = 123
  )

  cat("\nBinary organic/isotopic LDA test accuracy:\n")
  print(binary_results$test_accuracy)
} else {
  cat("\nBinary LDA skipped because the binary grouping could not be created.\n")
}


# ------------------------------------------------------------
# 13. Repeated train/test validation
# ------------------------------------------------------------

run_single_split_accuracy <- function(data,
                                      outcome_col,
                                      predictor_vars,
                                      p_train = 0.60,
                                      seed = 123,
                                      permute_labels = FALSE) {
  set.seed(seed)

  split_data <- data %>%
    dplyr::select(
      dplyr::all_of(outcome_col),
      dplyr::all_of(predictor_vars)
    ) %>%
    tidyr::drop_na()

  colnames(split_data)[1] <- "Outcome"
  split_data$Outcome <- droplevels(as.factor(split_data$Outcome))

  if (permute_labels) {
    split_data$Outcome <- sample(split_data$Outcome)
    split_data$Outcome <- droplevels(as.factor(split_data$Outcome))
  }

  split_index <- caret::createDataPartition(
    split_data$Outcome,
    p = p_train,
    list = FALSE
  )

  training_data <- split_data[split_index, , drop = FALSE]
  test_data <- split_data[-split_index, , drop = FALSE]

  y_train <- droplevels(as.factor(training_data$Outcome))
  y_test <- factor(test_data$Outcome, levels = levels(y_train))

  X_train_raw <- training_data %>%
    dplyr::select(dplyr::all_of(predictor_vars))

  X_test_raw <- test_data %>%
    dplyr::select(dplyr::all_of(predictor_vars))

  scaled_data <- autoscale_train_test(X_train_raw, X_test_raw)

  X_train <- scaled_data$train
  X_test <- scaled_data$test

  original_names <- colnames(X_train)
  model_names <- make.names(original_names, unique = TRUE)

  colnames(X_train) <- model_names
  colnames(X_test) <- model_names

  fit <- tryCatch(
    {
      lda_model <- MASS::lda(
        x = X_train,
        grouping = y_train
      )

      prediction <- predict(lda_model, newdata = X_test)
      predicted_class <- factor(prediction$class, levels = levels(y_train))
      accuracy <- mean(predicted_class == y_test)

      coefficients <- as.data.frame(lda_model$scaling)
      coefficients$Model_variable <- rownames(coefficients)

      variable_map <- data.frame(
        Model_variable = model_names,
        Original_variable = original_names
      )

      coefficients <- coefficients %>%
        dplyr::left_join(variable_map, by = "Model_variable")

      ld_columns <- grep("^LD", names(coefficients), value = TRUE)

      coefficients$Global_importance <- rowSums(
        abs(coefficients[, ld_columns, drop = FALSE])
      )

      coefficients <- coefficients %>%
        dplyr::arrange(dplyr::desc(Global_importance)) %>%
        dplyr::mutate(
          Seed = seed,
          Rank = dplyr::row_number()
        )

      list(
        success = TRUE,
        metrics = data.frame(
          Seed = seed,
          Accuracy = accuracy,
          Training_n = nrow(training_data),
          Test_n = nrow(test_data),
          Permuted = permute_labels
        ),
        importance = coefficients
      )
    },
    error = function(e) {
      message("Split failed for seed ", seed, ": ", e$message)

      list(
        success = FALSE,
        metrics = data.frame(
          Seed = seed,
          Accuracy = NA_real_,
          Training_n = NA_integer_,
          Test_n = NA_integer_,
          Permuted = permute_labels
        ),
        importance = NULL
      )
    }
  )

  fit
}

run_repeated_validation <- function(data,
                                    outcome_col,
                                    predictor_vars,
                                    model_name,
                                    output_prefix,
                                    n_splits = 100,
                                    n_perm = 500,
                                    p_train = 0.60) {
  cat("\nRunning repeated train/test validation for ", model_name, "...\n", sep = "")

  repeated_results <- vector("list", n_splits)

  for (i in seq_len(n_splits)) {
    repeated_results[[i]] <- run_single_split_accuracy(
      data = data,
      outcome_col = outcome_col,
      predictor_vars = predictor_vars,
      p_train = p_train,
      seed = 1000 + i,
      permute_labels = FALSE
    )
  }

  repeated_metrics <- purrr::map_dfr(repeated_results, "metrics")
  repeated_importance <- purrr::map_dfr(repeated_results, "importance")

  repeated_summary <- repeated_metrics %>%
    dplyr::filter(!is.na(Accuracy)) %>%
    dplyr::summarise(
      Valid_splits = dplyr::n(),
      Mean_accuracy = mean(Accuracy),
      SD_accuracy = sd(Accuracy),
      Minimum_accuracy = min(Accuracy),
      Maximum_accuracy = max(Accuracy),
      Mean_accuracy_percent = 100 * Mean_accuracy,
      SD_accuracy_percent = 100 * SD_accuracy,
      Minimum_accuracy_percent = 100 * Minimum_accuracy,
      Maximum_accuracy_percent = 100 * Maximum_accuracy
    )

  write.csv(
    repeated_metrics,
    file.path(tables_dir, paste0(output_prefix, "_repeated_split_metrics.csv")),
    row.names = FALSE
  )

  write.csv(
    repeated_summary,
    file.path(tables_dir, paste0(output_prefix, "_repeated_split_summary.csv")),
    row.names = FALSE
  )

  repeated_plot <- repeated_metrics %>%
    dplyr::filter(!is.na(Accuracy)) %>%
    ggplot2::ggplot(ggplot2::aes(x = model_name, y = Accuracy * 100)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.5) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::labs(
      title = paste0(model_name, " - repeated train/test validation"),
      x = NULL,
      y = "Accuracy (%)"
    )

  ggplot2::ggsave(
    file.path(figures_dir, paste0(output_prefix, "_repeated_split_accuracy.png")),
    repeated_plot,
    width = 6,
    height = 5,
    dpi = 300
  )

  observed_mean_accuracy <- repeated_summary$Mean_accuracy[1]

  cat("\nRunning permutation test for ", model_name, "...\n", sep = "")

  permutation_results <- vector("list", n_perm)

  for (i in seq_len(n_perm)) {
    permutation_results[[i]] <- run_single_split_accuracy(
      data = data,
      outcome_col = outcome_col,
      predictor_vars = predictor_vars,
      p_train = p_train,
      seed = 5000 + i,
      permute_labels = TRUE
    )
  }

  permutation_metrics <- purrr::map_dfr(permutation_results, "metrics") %>%
    dplyr::filter(!is.na(Accuracy))

  empirical_p_value <- (
    sum(permutation_metrics$Accuracy >= observed_mean_accuracy) + 1
  ) / (nrow(permutation_metrics) + 1)

  permutation_summary <- data.frame(
    Observed_mean_accuracy = observed_mean_accuracy,
    Permuted_mean_accuracy = mean(permutation_metrics$Accuracy),
    Permuted_SD_accuracy = sd(permutation_metrics$Accuracy),
    Permuted_minimum_accuracy = min(permutation_metrics$Accuracy),
    Permuted_maximum_accuracy = max(permutation_metrics$Accuracy),
    Empirical_p_value = empirical_p_value
  )

  write.csv(
    permutation_metrics,
    file.path(tables_dir, paste0(output_prefix, "_permutation_metrics.csv")),
    row.names = FALSE
  )

  write.csv(
    permutation_summary,
    file.path(tables_dir, paste0(output_prefix, "_permutation_summary.csv")),
    row.names = FALSE
  )

  permutation_plot <- permutation_metrics %>%
    ggplot2::ggplot(ggplot2::aes(x = Accuracy * 100)) +
    ggplot2::geom_histogram(bins = 30, color = "black") +
    ggplot2::geom_vline(
      xintercept = observed_mean_accuracy * 100,
      linetype = "dashed",
      linewidth = 1
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::labs(
      title = paste0(model_name, " - permutation test"),
      x = "Accuracy with permuted labels (%)",
      y = "Frequency"
    )

  ggplot2::ggsave(
    file.path(figures_dir, paste0(output_prefix, "_permutation_test.png")),
    permutation_plot,
    width = 7,
    height = 5,
    dpi = 300
  )

  variable_stability <- NULL

  if (nrow(repeated_importance) > 0) {
    top_n <- min(10, length(predictor_vars))

    variable_stability <- repeated_importance %>%
      dplyr::group_by(Seed) %>%
      dplyr::arrange(Rank, .by_group = TRUE) %>%
      dplyr::slice_head(n = top_n) %>%
      dplyr::ungroup() %>%
      dplyr::count(Original_variable, name = "Top_n_frequency") %>%
      dplyr::mutate(
        Top_n_frequency_percent = 100 * Top_n_frequency / n_splits
      ) %>%
      dplyr::arrange(dplyr::desc(Top_n_frequency_percent))

    write.csv(
      variable_stability,
      file.path(tables_dir, paste0(output_prefix, "_variable_importance_stability_top10.csv")),
      row.names = FALSE
    )

    stability_plot <- variable_stability %>%
      ggplot2::ggplot(
        ggplot2::aes(
          x = reorder(Original_variable, Top_n_frequency_percent),
          y = Top_n_frequency_percent
        )
      ) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::theme_bw(base_size = 14) +
      ggplot2::labs(
        title = paste0(model_name, " - variable-importance stability"),
        x = "Variable",
        y = paste0("Frequency in top ", top_n, " variables (%)")
      )

    ggplot2::ggsave(
      file.path(figures_dir, paste0(output_prefix, "_variable_importance_stability_top10.png")),
      stability_plot,
      width = 8,
      height = 6,
      dpi = 300
    )
  }

  list(
    repeated_metrics = repeated_metrics,
    repeated_summary = repeated_summary,
    permutation_metrics = permutation_metrics,
    permutation_summary = permutation_summary,
    variable_stability = variable_stability
  )
}


# ------------------------------------------------------------
# 14. Robustness evaluation
# ------------------------------------------------------------

if (nlevels(clean_data$Group) > 2) {
  multiclass_robustness <- run_repeated_validation(
    data = clean_data,
    outcome_col = "Group",
    predictor_vars = organic_isotopic_vars,
    model_name = "Organic/isotopic multiclass LDA",
    output_prefix = "organic_isotopic_multiclass_lda",
    n_splits = n_repeated_splits,
    n_perm = n_permutations,
    p_train = training_fraction
  )
}

if (nlevels(droplevels(binary_data$Binary_Group)) == 2) {
  binary_robustness <- run_repeated_validation(
    data = binary_data,
    outcome_col = "Binary_Group",
    predictor_vars = organic_isotopic_vars,
    model_name = "Organic/isotopic binary LDA",
    output_prefix = "organic_isotopic_binary_lda",
    n_splits = n_repeated_splits,
    n_perm = n_permutations,
    p_train = training_fraction
  )
}


# ------------------------------------------------------------
# 15. Session information
# ------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo_organic_isotopic_model.txt")
)

cat("\nOrganic/isotopic LDA workflow completed successfully.\n")
cat("Tables were saved in: ", tables_dir, "\n", sep = "")
cat("Figures were saved in: ", figures_dir, "\n", sep = "")
