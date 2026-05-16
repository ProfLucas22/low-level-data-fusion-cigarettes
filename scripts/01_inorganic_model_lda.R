# ============================================================
# INORGANIC LDA MODEL FOR COMMERCIAL CIGARETTE CLASSIFICATION
# Manuscript supplementary reference material
# ============================================================
#
# Description:
# This script performs the chemometric workflow for the inorganic descriptor
# block used in the classification of commercial cigarette samples.
#
# The workflow evaluates two classification tasks:
#   1. Multiclass classification
#   2. Binary classification: Multinational vs Others
#
# Main steps:
#   1. Load packages and input data
#   2. Assign or standardize class labels
#   3. Select inorganic descriptors
#   4. Generate descriptive statistics and plots
#   5. Apply class-wise ROBPCA-based multivariate outlier screening
#   6. Split the cleaned dataset into training and test sets
#   7. Apply log1p transformation and autoscaling using training-set
#      parameters only
#   8. Fit LDA models
#   9. Generate confusion matrices, LD score plots, ROC/AUC results,
#      variable-importance tables, and robustness summaries
#
# Expected input:
#   Option 1: an object named teste_priliminar_tese_2 already loaded in R
#   Option 2: a CSV file named inorganic_data.csv in data/processed/ or data/
#   Option 3: a CSV file named teste_priliminar_tese_2.csv in data/processed/ or data/
#
# Required inorganic variables:
#   Zn, Sr, Ni, Mn, Fe, Cu, Co, Cd, Ba, B
#
# If the input file does not contain a column named "grupos", class labels are
# assigned using the sequence originally used in the manuscript workflow.
#
# ============================================================


# ------------------------------------------------------------
# 1. Global settings
# ------------------------------------------------------------

set.seed(123)

RUN_ROBUSTNESS <- TRUE
N_REPETITIONS <- 100
N_PERMUTATIONS <- 500
TRAINING_PROPORTION <- 0.60

output_dir <- file.path("outputs", "inorganic_model")
tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)


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
# 3. Helper functions
# ------------------------------------------------------------

translate_group_labels <- function(x) {
  x <- as.character(x)
  dplyr::recode(
    x,
    "Multinacional" = "Multinational",
    "Contrabandeado" = "Smuggled",
    "Falsificado" = "Counterfeit",
    "Sem_Registro_Sanitário" = "No_Health_Registration",
    "Sem_Registro_Sanitario" = "No_Health_Registration",
    "Regional" = "Regional",
    "Outros" = "Others",
    .default = x
  )
}

assign_default_multiclass_labels <- function(n) {
  labels <- c(
    rep("Regional", 3),
    rep("Falsificado", 3),
    rep("Contrabandeado", 3),
    rep("Regional", 2),
    rep("Sem_Registro_Sanitário", 3),
    rep("Contrabandeado", 6),
    rep("Falsificado", 9),
    rep("Contrabandeado", 3),
    rep("Sem_Registro_Sanitário", 3),
    rep("Contrabandeado", 3),
    rep("Regional", 3),
    rep("Multinacional", 3),
    rep("Regional", 6),
    rep("Sem_Registro_Sanitário", 3),
    rep("Regional", 5),
    rep("Contrabandeado", 8),
    rep("Multinacional", 25),
    "Regional",
    rep("Multinacional", 2),
    rep("Falsificado", 3),
    rep("Regional", 3),
    rep("Contrabandeado", 3),
    rep("Regional", 4)
  )

  if (length(labels) != n) {
    stop(
      "The default multiclass label vector has length ", length(labels),
      ", but the dataset has ", n, " rows. Please provide a 'grupos' column."
    )
  }

  translate_group_labels(labels)
}

create_binary_labels <- function(multiclass_labels) {
  ifelse(multiclass_labels == "Multinational", "Multinational", "Others")
}

safe_log1p <- function(x) {
  x <- as.numeric(x)
  if (any(x < -1, na.rm = TRUE)) {
    stop("log1p transformation cannot be applied because values lower than -1 were found.")
  }
  log1p(x)
}

remove_zero_variance_vars <- function(data, vars) {
  vars[sapply(data[vars], function(x) {
    x <- as.numeric(x)
    x <- x[!is.na(x)]
    length(unique(x)) > 1 && stats::sd(x) > 0
  })]
}

write_table <- function(x, filename) {
  utils::write.csv2(
    x,
    file.path(tables_dir, filename),
    row.names = FALSE
  )
}

save_ggplot <- function(plot_object, filename, width = 8, height = 5) {
  ggplot2::ggsave(
    filename = file.path(figures_dir, filename),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300
  )
}

prepare_train_test_matrices <- function(training_data,
                                        test_data,
                                        predictor_vars,
                                        group_col = "group") {
  X_train_raw <- training_data %>%
    dplyr::select(dplyr::all_of(predictor_vars)) %>%
    as.data.frame()

  X_test_raw <- test_data %>%
    dplyr::select(dplyr::all_of(predictor_vars)) %>%
    as.data.frame()

  X_train_log <- as.data.frame(lapply(X_train_raw, safe_log1p))
  X_test_log <- as.data.frame(lapply(X_test_raw, safe_log1p))

  colnames(X_train_log) <- predictor_vars
  colnames(X_test_log) <- predictor_vars

  vars_ok <- remove_zero_variance_vars(X_train_log, predictor_vars)

  if (length(vars_ok) < 2) {
    stop("Fewer than two predictors with variability were retained in the training set.")
  }

  X_train_log <- X_train_log[, vars_ok, drop = FALSE]
  X_test_log <- X_test_log[, vars_ok, drop = FALSE]

  centers <- sapply(X_train_log, mean, na.rm = TRUE)
  scales <- sapply(X_train_log, stats::sd, na.rm = TRUE)
  scales[is.na(scales) | scales == 0] <- 1

  X_train_scaled <- sweep(X_train_log, 2, centers, "-")
  X_train_scaled <- sweep(X_train_scaled, 2, scales, "/")

  X_test_scaled <- sweep(X_test_log, 2, centers, "-")
  X_test_scaled <- sweep(X_test_scaled, 2, scales, "/")

  model_names <- make.names(colnames(X_train_scaled), unique = TRUE)
  variable_map <- data.frame(
    model_variable = model_names,
    original_variable = colnames(X_train_scaled),
    stringsAsFactors = FALSE
  )

  colnames(X_train_scaled) <- model_names
  colnames(X_test_scaled) <- model_names

  list(
    X_train = as.data.frame(X_train_scaled),
    X_test = as.data.frame(X_test_scaled),
    selected_vars = vars_ok,
    variable_map = variable_map,
    centers = centers,
    scales = scales
  )
}

detect_classwise_robpca_outliers <- function(data,
                                             predictor_vars,
                                             group_col = "group",
                                             k_max = 2) {
  data_for_robpca <- data %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(predictor_vars),
        safe_log1p
      )
    )

  detect_one_group <- function(df, vars) {
    out <- df

    out$robpca_flag <- NA_integer_
    out$robpca_sd <- NA_real_
    out$robpca_od <- NA_real_
    out$robpca_cutoff_sd <- NA_real_
    out$robpca_cutoff_od <- NA_real_
    out$robpca_outlier <- NA
    out$robpca_note <- NA_character_

    X <- out[, vars, drop = FALSE] %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric))

    vars_keep <- remove_zero_variance_vars(X, vars)
    X <- X[, vars_keep, drop = FALSE]

    if (ncol(X) < 2) {
      out$robpca_note <- "Fewer than two variables with variability in this class."
      return(out)
    }

    complete_rows <- complete.cases(X)

    if (sum(complete_rows) < 4) {
      out$robpca_note <- "Fewer than four complete observations in this class."
      return(out)
    }

    X_complete <- as.matrix(X[complete_rows, , drop = FALSE])
    X_scaled <- scale(X_complete)

    n_components <- min(k_max, ncol(X_scaled) - 1)

    if (n_components < 1) {
      out$robpca_note <- "Insufficient dimensionality for ROBPCA."
      return(out)
    }

    fit <- tryCatch(
      rrcov::PcaHubert(X_scaled, k = n_components, scale = FALSE),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      out$robpca_note <- paste("ROBPCA error:", fit$message)
      return(out)
    }

    out$robpca_flag[complete_rows] <- as.integer(fit@flag)
    out$robpca_sd[complete_rows] <- as.numeric(fit@sd)
    out$robpca_od[complete_rows] <- as.numeric(fit@od)
    out$robpca_cutoff_sd[complete_rows] <- as.numeric(fit@cutoff.sd)
    out$robpca_cutoff_od[complete_rows] <- as.numeric(fit@cutoff.od)
    out$robpca_outlier[complete_rows] <- out$robpca_flag[complete_rows] == 0
    out$robpca_note <- paste0("ROBPCA applied using ", ncol(X_scaled), " variables.")

    out
  }

  data_for_robpca %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::group_modify(~ detect_one_group(.x, predictor_vars)) %>%
    dplyr::ungroup()
}

compute_class_metrics <- function(confusion_table) {
  classes <- rownames(confusion_table)

  metrics <- data.frame(
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

run_lda_workflow <- function(data,
                             group_col,
                             predictor_vars,
                             model_label,
                             output_prefix,
                             seed = 123,
                             training_proportion = 0.60) {
  set.seed(seed)

  model_data <- data %>%
    dplyr::select(
      dplyr::all_of(group_col),
      dplyr::all_of(predictor_vars)
    ) %>%
    tidyr::drop_na()

  colnames(model_data)[colnames(model_data) == group_col] <- "group"
  model_data$group <- as.factor(model_data$group)

  class_counts <- table(model_data$group)
  if (any(class_counts < 2)) {
    stop("At least one class has fewer than two observations. LDA cannot be fitted reliably.")
  }

  idx_train <- caret::createDataPartition(
    model_data$group,
    p = training_proportion,
    list = FALSE
  )

  training_data <- model_data[idx_train, , drop = FALSE]
  test_data <- model_data[-idx_train, , drop = FALSE]

  y_train <- as.factor(training_data$group)
  y_test <- factor(test_data$group, levels = levels(y_train))

  preprocessing <- prepare_train_test_matrices(
    training_data = training_data,
    test_data = test_data,
    predictor_vars = predictor_vars,
    group_col = "group"
  )

  X_train <- preprocessing$X_train
  X_test <- preprocessing$X_test

  lda_model <- MASS::lda(
    x = X_train,
    grouping = y_train
  )

  pred_train <- predict(lda_model, newdata = X_train)
  pred_test <- predict(lda_model, newdata = X_test)

  class_train <- factor(pred_train$class, levels = levels(y_train))
  class_test <- factor(pred_test$class, levels = levels(y_train))

  posterior_train <- as.data.frame(pred_train$posterior)
  posterior_test <- as.data.frame(pred_test$posterior)

  scores_train <- as.data.frame(pred_train$x)
  scores_test <- as.data.frame(pred_test$x)

  cm_train <- table(Observed = y_train, Predicted = class_train)
  cm_test <- table(Observed = y_test, Predicted = class_test)

  acc_train <- mean(class_train == y_train)
  acc_test <- mean(class_test == y_test)

  train_metrics <- compute_class_metrics(cm_train) %>%
    dplyr::mutate(Set = "Training")

  test_metrics <- compute_class_metrics(cm_test) %>%
    dplyr::mutate(Set = "Test")

  all_metrics <- dplyr::bind_rows(train_metrics, test_metrics)

  overall_metrics <- data.frame(
    Model = model_label,
    Training_accuracy = acc_train,
    Test_accuracy = acc_test,
    N_training = nrow(training_data),
    N_test = nrow(test_data)
  )

  # Variable importance
  coef_lda <- as.data.frame(lda_model$scaling)
  coef_lda$model_variable <- rownames(coef_lda)

  coef_lda <- coef_lda %>%
    dplyr::left_join(preprocessing$variable_map, by = "model_variable")

  ld_cols <- grep("^LD", names(coef_lda), value = TRUE)

  coef_lda$Global_importance <- rowSums(
    abs(coef_lda[, ld_cols, drop = FALSE])
  )

  coef_lda <- coef_lda %>%
    dplyr::arrange(dplyr::desc(Global_importance))

  # Explained variance
  eigenvalues <- lda_model$svd^2
  explained_variance <- eigenvalues / sum(eigenvalues) * 100

  explained_variance_df <- data.frame(
    LD = paste0("LD", seq_along(explained_variance)),
    Eigenvalue = eigenvalues,
    Explained_variance_percent = explained_variance,
    Cumulative_variance_percent = cumsum(explained_variance)
  )

  # Scores
  scores_train$group <- y_train
  scores_train$Predicted <- class_train
  scores_train$Set <- "Training"

  scores_test$group <- y_test
  scores_test$Predicted <- class_test
  scores_test$Set <- "Test"

  scores_all <- dplyr::bind_rows(scores_train, scores_test)

  # ROC/AUC one-vs-rest
  classes <- levels(y_train)
  auc_values <- numeric(length(classes))
  names(auc_values) <- classes
  roc_list <- list()

  for (class_i in classes) {
    if (!class_i %in% colnames(posterior_test)) {
      auc_values[class_i] <- NA_real_
      next
    }

    y_binary <- ifelse(y_test == class_i, 1, 0)

    if (length(unique(y_binary)) < 2) {
      auc_values[class_i] <- NA_real_
      next
    }

    roc_i <- pROC::roc(
      response = y_binary,
      predictor = posterior_test[[class_i]],
      quiet = TRUE
    )

    roc_list[[class_i]] <- roc_i
    auc_values[class_i] <- as.numeric(pROC::auc(roc_i))
  }

  auc_df <- data.frame(
    Class = names(auc_values),
    AUC_one_vs_rest = as.numeric(auc_values)
  )

  # Export tables
  write_table(as.data.frame(cm_train), paste0(output_prefix, "_confusion_matrix_training.csv"))
  write_table(as.data.frame(cm_test), paste0(output_prefix, "_confusion_matrix_test.csv"))
  write_table(overall_metrics, paste0(output_prefix, "_overall_metrics.csv"))
  write_table(all_metrics, paste0(output_prefix, "_class_metrics.csv"))
  write_table(coef_lda, paste0(output_prefix, "_variable_importance.csv"))
  write_table(explained_variance_df, paste0(output_prefix, "_explained_variance.csv"))
  write_table(auc_df, paste0(output_prefix, "_auc_one_vs_rest.csv"))
  write_table(scores_all, paste0(output_prefix, "_lda_scores.csv"))

  # Confusion matrix plot - test set
  cm_test_df <- as.data.frame(cm_test)

  p_cm <- ggplot2::ggplot(
    cm_test_df,
    ggplot2::aes(x = Predicted, y = Observed, fill = Freq)
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = Freq), size = 5) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::labs(
      title = paste0("Test confusion matrix - ", model_label),
      x = "Predicted class",
      y = "Observed class",
      fill = "N"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  save_ggplot(p_cm, paste0(output_prefix, "_confusion_matrix_test.png"), width = 7, height = 5)

  # LD score plot
  if ("LD1" %in% names(scores_all) && "LD2" %in% names(scores_all)) {
    x_label <- paste0("LD1 (", round(explained_variance[1], 2), "%)")
    y_label <- paste0("LD2 (", round(explained_variance[2], 2), "%)")

    p_scores <- ggplot2::ggplot(
      scores_all,
      ggplot2::aes(x = LD1, y = LD2, color = group, shape = Set)
    ) +
      ggplot2::geom_point(size = 3, alpha = 0.85) +
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::labs(
        title = paste0("LDA score plot - ", model_label),
        x = x_label,
        y = y_label,
        color = "Class",
        shape = "Set"
      )
  } else if ("LD1" %in% names(scores_all)) {
    x_label <- paste0("LD1 (", round(explained_variance[1], 2), "%)")

    p_scores <- ggplot2::ggplot(
      scores_all,
      ggplot2::aes(x = LD1, y = 0, color = group, shape = Set)
    ) +
      ggplot2::geom_point(size = 3, alpha = 0.85) +
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::labs(
        title = paste0("LDA score plot - ", model_label),
        x = x_label,
        y = "",
        color = "Class",
        shape = "Set"
      ) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
  } else {
    p_scores <- NULL
  }

  if (!is.null(p_scores)) {
    save_ggplot(p_scores, paste0(output_prefix, "_lda_scores.png"), width = 8, height = 5)
  }

  # Explained variance plot
  p_var <- ggplot2::ggplot(
    explained_variance_df,
    ggplot2::aes(x = LD, y = Explained_variance_percent, group = 1)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(
      ggplot2::aes(label = round(Explained_variance_percent, 2)),
      vjust = -0.8,
      size = 3.5
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::labs(
      title = paste0("Explained variance - ", model_label),
      x = "Linear discriminant",
      y = "Explained variance (%)"
    )

  save_ggplot(p_var, paste0(output_prefix, "_explained_variance.png"), width = 7, height = 5)

  # Variable importance plot
  top_variables <- coef_lda %>%
    dplyr::slice_head(n = min(10, nrow(coef_lda)))

  p_importance <- ggplot2::ggplot(
    top_variables,
    ggplot2::aes(
      x = reorder(original_variable, Global_importance),
      y = Global_importance
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::labs(
      title = paste0("Variable importance - ", model_label),
      x = "Variable",
      y = "Global importance"
    )

  save_ggplot(p_importance, paste0(output_prefix, "_variable_importance.png"), width = 7, height = 5)

  # ROC plot
  if (length(roc_list) > 0) {
    png(
      filename = file.path(figures_dir, paste0(output_prefix, "_roc_one_vs_rest.png")),
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
      main = paste0("One-vs-rest ROC curves - ", model_label)
    )

    abline(a = 1, b = -1, lty = 2, col = "gray50")

    colors <- grDevices::rainbow(length(roc_list))

    for (i in seq_along(roc_list)) {
      plot(
        roc_list[[i]],
        add = TRUE,
        col = colors[i],
        lwd = 2,
        legacy.axes = FALSE
      )
    }

    legend(
      "bottomleft",
      legend = paste0(
        names(roc_list),
        " (AUC = ",
        round(auc_df$AUC_one_vs_rest[match(names(roc_list), auc_df$Class)], 3),
        ")"
      ),
      col = colors,
      lwd = 2,
      cex = 0.8,
      bty = "n"
    )

    dev.off()
  }

  list(
    model = lda_model,
    preprocessing = preprocessing,
    training_data = training_data,
    test_data = test_data,
    predictions_train = pred_train,
    predictions_test = pred_test,
    confusion_training = cm_train,
    confusion_test = cm_test,
    overall_metrics = overall_metrics,
    class_metrics = all_metrics,
    variable_importance = coef_lda,
    explained_variance = explained_variance_df,
    auc = auc_df,
    scores = scores_all
  )
}

run_single_split_accuracy <- function(data,
                                      group_col,
                                      predictor_vars,
                                      seed,
                                      training_proportion = 0.60,
                                      permute_labels = FALSE) {
  set.seed(seed)

  model_data <- data %>%
    dplyr::select(
      dplyr::all_of(group_col),
      dplyr::all_of(predictor_vars)
    ) %>%
    tidyr::drop_na()

  colnames(model_data)[colnames(model_data) == group_col] <- "group"
  model_data$group <- as.factor(model_data$group)

  if (permute_labels) {
    model_data$group <- sample(model_data$group)
    model_data$group <- as.factor(model_data$group)
  }

  result <- tryCatch({
    idx_train <- caret::createDataPartition(
      model_data$group,
      p = training_proportion,
      list = FALSE
    )

    training_data <- model_data[idx_train, , drop = FALSE]
    test_data <- model_data[-idx_train, , drop = FALSE]

    y_train <- as.factor(training_data$group)
    y_test <- factor(test_data$group, levels = levels(y_train))

    preprocessing <- prepare_train_test_matrices(
      training_data = training_data,
      test_data = test_data,
      predictor_vars = predictor_vars,
      group_col = "group"
    )

    lda_model <- MASS::lda(
      x = preprocessing$X_train,
      grouping = y_train
    )

    pred <- predict(lda_model, newdata = preprocessing$X_test)
    pred_class <- factor(pred$class, levels = levels(y_train))

    accuracy <- mean(pred_class == y_test)

    coef_lda <- as.data.frame(lda_model$scaling)
    coef_lda$model_variable <- rownames(coef_lda)

    coef_lda <- coef_lda %>%
      dplyr::left_join(preprocessing$variable_map, by = "model_variable")

    ld_cols <- grep("^LD", names(coef_lda), value = TRUE)

    coef_lda$Global_importance <- rowSums(
      abs(coef_lda[, ld_cols, drop = FALSE])
    )

    coef_lda <- coef_lda %>%
      dplyr::arrange(dplyr::desc(Global_importance)) %>%
      dplyr::mutate(
        seed = seed,
        rank = dplyr::row_number()
      )

    list(
      metrics = data.frame(
        seed = seed,
        accuracy = accuracy,
        n_training = nrow(training_data),
        n_test = nrow(test_data),
        permuted = permute_labels
      ),
      variable_importance = coef_lda
    )
  }, error = function(e) {
    message("Split ", seed, " failed: ", e$message)
    list(
      metrics = data.frame(
        seed = seed,
        accuracy = NA_real_,
        n_training = NA_integer_,
        n_test = NA_integer_,
        permuted = permute_labels
      ),
      variable_importance = NULL
    )
  })

  result
}

run_robustness_analysis <- function(data,
                                    group_col,
                                    predictor_vars,
                                    output_prefix,
                                    n_repetitions = 100,
                                    n_permutations = 500,
                                    training_proportion = 0.60) {
  repeated_results <- vector("list", n_repetitions)

  for (i in seq_len(n_repetitions)) {
    repeated_results[[i]] <- run_single_split_accuracy(
      data = data,
      group_col = group_col,
      predictor_vars = predictor_vars,
      seed = 1000 + i,
      training_proportion = training_proportion,
      permute_labels = FALSE
    )
  }

  repeated_metrics <- purrr::map_dfr(repeated_results, "metrics") %>%
    dplyr::filter(!is.na(accuracy))

  repeated_importance <- purrr::map_dfr(repeated_results, "variable_importance")

  robustness_summary <- repeated_metrics %>%
    dplyr::summarise(
      valid_repetitions = dplyr::n(),
      mean_accuracy = mean(accuracy),
      sd_accuracy = stats::sd(accuracy),
      min_accuracy = min(accuracy),
      max_accuracy = max(accuracy)
    ) %>%
    dplyr::mutate(
      mean_accuracy_percent = 100 * mean_accuracy,
      sd_accuracy_percent = 100 * sd_accuracy,
      min_accuracy_percent = 100 * min_accuracy,
      max_accuracy_percent = 100 * max_accuracy
    )

  write_table(repeated_metrics, paste0(output_prefix, "_robustness_repeated_splits.csv"))
  write_table(robustness_summary, paste0(output_prefix, "_robustness_summary.csv"))

  p_repeated <- ggplot2::ggplot(
    repeated_metrics,
    ggplot2::aes(x = "Repeated splits", y = accuracy * 100)
  ) +
    ggplot2::geom_boxplot(width = 0.35) +
    ggplot2::geom_jitter(width = 0.08, alpha = 0.5) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::labs(
      title = "Repeated train/test validation",
      x = "",
      y = "Accuracy (%)"
    )

  save_ggplot(
    p_repeated,
    paste0(output_prefix, "_robustness_repeated_splits.png"),
    width = 6,
    height = 5
  )

  if (n_permutations > 0) {
    permutation_results <- vector("list", n_permutations)

    for (i in seq_len(n_permutations)) {
      permutation_results[[i]] <- run_single_split_accuracy(
        data = data,
        group_col = group_col,
        predictor_vars = predictor_vars,
        seed = 5000 + i,
        training_proportion = training_proportion,
        permute_labels = TRUE
      )
    }

    permutation_metrics <- purrr::map_dfr(permutation_results, "metrics") %>%
      dplyr::filter(!is.na(accuracy))

    observed_accuracy <- robustness_summary$mean_accuracy

    empirical_p <- (sum(permutation_metrics$accuracy >= observed_accuracy) + 1) /
      (nrow(permutation_metrics) + 1)

    permutation_summary <- data.frame(
      observed_mean_accuracy = observed_accuracy,
      permuted_mean_accuracy = mean(permutation_metrics$accuracy),
      permuted_sd_accuracy = stats::sd(permutation_metrics$accuracy),
      permuted_min_accuracy = min(permutation_metrics$accuracy),
      permuted_max_accuracy = max(permutation_metrics$accuracy),
      empirical_p_value = empirical_p
    )

    write_table(permutation_metrics, paste0(output_prefix, "_permutation_metrics.csv"))
    write_table(permutation_summary, paste0(output_prefix, "_permutation_summary.csv"))

    p_perm <- ggplot2::ggplot(
      permutation_metrics,
      ggplot2::aes(x = accuracy * 100)
    ) +
      ggplot2::geom_histogram(bins = 30, color = "black") +
      ggplot2::geom_vline(
        xintercept = observed_accuracy * 100,
        linetype = "dashed",
        linewidth = 1
      ) +
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::labs(
        title = "Permutation test",
        x = "Accuracy with permuted labels (%)",
        y = "Frequency"
      )

    save_ggplot(
      p_perm,
      paste0(output_prefix, "_permutation_test.png"),
      width = 7,
      height = 5
    )
  }

  if (nrow(repeated_importance) > 0) {
    top_n <- min(10, length(predictor_vars))

    variable_stability <- repeated_importance %>%
      dplyr::group_by(seed) %>%
      dplyr::arrange(rank, .by_group = TRUE) %>%
      dplyr::slice_head(n = top_n) %>%
      dplyr::ungroup() %>%
      dplyr::count(original_variable, name = "frequency_top_n") %>%
      dplyr::mutate(
        frequency_percent = 100 * frequency_top_n / n_repetitions
      ) %>%
      dplyr::arrange(dplyr::desc(frequency_percent))

    write_table(variable_stability, paste0(output_prefix, "_variable_importance_stability.csv"))

    p_stability <- ggplot2::ggplot(
      variable_stability,
      ggplot2::aes(
        x = reorder(original_variable, frequency_percent),
        y = frequency_percent
      )
    ) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::labs(
        title = paste0("Variable-importance stability - top ", top_n),
        x = "Variable",
        y = paste0("Frequency in top ", top_n, " variables (%)")
      )

    save_ggplot(
      p_stability,
      paste0(output_prefix, "_variable_importance_stability.png"),
      width = 7,
      height = 5
    )
  }

  invisible(TRUE)
}


# ------------------------------------------------------------
# 4. Load data
# ------------------------------------------------------------

if (exists("teste_priliminar_tese_2")) {
  raw_data <- teste_priliminar_tese_2
} else if (file.exists("data/processed/inorganic_data.csv")) {
  raw_data <- read.csv(
    "data/processed/inorganic_data.csv",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
} else if (file.exists("data/inorganic_data.csv")) {
  raw_data <- read.csv(
    "data/inorganic_data.csv",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
} else if (file.exists("data/processed/teste_priliminar_tese_2.csv")) {
  raw_data <- read.csv(
    "data/processed/teste_priliminar_tese_2.csv",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
} else if (file.exists("data/teste_priliminar_tese_2.csv")) {
  raw_data <- read.csv(
    "data/teste_priliminar_tese_2.csv",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
} else {
  stop(
    "Input data not found. Please load teste_priliminar_tese_2 in R, ",
    "or place inorganic_data.csv or teste_priliminar_tese_2.csv in data/processed/ or data/."
  )
}


# ------------------------------------------------------------
# 5. Define variables and class labels
# ------------------------------------------------------------

inorganic_vars <- c("Zn", "Sr", "Ni", "Mn", "Fe", "Cu", "Co", "Cd", "Ba", "B")

missing_vars <- setdiff(inorganic_vars, names(raw_data))

if (length(missing_vars) > 0) {
  stop(
    "The following inorganic variables are missing from the dataset: ",
    paste(missing_vars, collapse = ", ")
  )
}

data_all <- raw_data %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(inorganic_vars),
      as.numeric
    )
  )

if ("grupos" %in% names(data_all)) {
  data_all$group_multiclass <- translate_group_labels(data_all$grupos)
} else if ("group" %in% names(data_all)) {
  data_all$group_multiclass <- translate_group_labels(data_all$group)
} else {
  data_all$group_multiclass <- assign_default_multiclass_labels(nrow(data_all))
}

data_all$group_multiclass <- factor(data_all$group_multiclass)

data_all$group_binary <- factor(
  create_binary_labels(data_all$group_multiclass),
  levels = c("Others", "Multinational")
)

cat("\nMulticlass class distribution in the original dataset:\n")
print(table(data_all$group_multiclass))

cat("\nBinary class distribution in the original dataset:\n")
print(table(data_all$group_binary))


# ------------------------------------------------------------
# 6. Descriptive statistics and plot
# ------------------------------------------------------------

descriptive_data <- data_all %>%
  dplyr::select(group_multiclass, dplyr::all_of(inorganic_vars)) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(inorganic_vars),
      safe_log1p
    )
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(inorganic_vars),
    names_to = "Variable",
    values_to = "Value"
  )

descriptive_summary <- descriptive_data %>%
  dplyr::group_by(group_multiclass, Variable) %>%
  dplyr::summarise(
    Mean = mean(Value, na.rm = TRUE),
    SD = stats::sd(Value, na.rm = TRUE),
    N = sum(!is.na(Value)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    SD = ifelse(is.na(SD), 0, SD),
    Label = paste0("(", round(Mean, 2), " ± ", round(SD, 2), ")")
  ) %>%
  dplyr::filter(!is.na(Mean))

write_table(descriptive_summary, "inorganic_descriptive_summary_log_scale.csv")

p_descriptive <- ggplot2::ggplot(
  descriptive_summary,
  ggplot2::aes(x = group_multiclass, y = Mean, color = group_multiclass)
) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = Mean - SD, ymax = Mean + SD),
    width = 0.12,
    linewidth = 0.5
  ) +
  ggrepel::geom_label_repel(
    ggplot2::aes(label = Label),
    size = 2.4,
    fill = "white",
    label.size = 0.15,
    label.padding = grid::unit(0.12, "lines"),
    box.padding = 0.25,
    point.padding = 0.15,
    force = 1.2,
    direction = "y",
    segment.color = "grey50",
    segment.size = 0.3,
    min.segment.length = 0,
    seed = 123,
    show.legend = FALSE
  ) +
  ggplot2::facet_wrap(~ Variable, scales = "free_y", ncol = 3) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0.05, 0.30))
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::labs(
    x = NULL,
    y = "Mean ± standard deviation (log1p scale)",
    color = "Class"
  ) +
  ggplot2::theme(
    legend.position = "bottom",
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = ggplot2::element_text(size = 9),
    strip.text = ggplot2::element_text(size = 10, face = "bold"),
    panel.spacing = grid::unit(1.3, "lines"),
    panel.grid.major = ggplot2::element_line(linewidth = 0.2),
    panel.grid.minor = ggplot2::element_blank()
  )

save_ggplot(
  p_descriptive,
  "inorganic_descriptive_summary_log_scale.png",
  width = 10,
  height = 7
)


# ------------------------------------------------------------
# 7. ROBPCA-based outlier screening
# ------------------------------------------------------------

multiclass_robpca <- detect_classwise_robpca_outliers(
  data = data_all %>%
    dplyr::select(group_multiclass, dplyr::all_of(inorganic_vars)),
  predictor_vars = inorganic_vars,
  group_col = "group_multiclass",
  k_max = 2
)

outlier_table_multiclass <- multiclass_robpca %>%
  dplyr::mutate(original_row = dplyr::row_number()) %>%
  dplyr::select(
    original_row,
    group_multiclass,
    robpca_flag,
    robpca_sd,
    robpca_od,
    robpca_cutoff_sd,
    robpca_cutoff_od,
    robpca_outlier,
    robpca_note
  )

write_table(outlier_table_multiclass, "inorganic_multiclass_robpca_outliers.csv")

cat("\nNumber of multiclass ROBPCA outliers:\n")
print(sum(outlier_table_multiclass$robpca_outlier, na.rm = TRUE))

cat("\nMulticlass ROBPCA outliers by class:\n")
print(table(
  Class = outlier_table_multiclass$group_multiclass,
  Outlier = outlier_table_multiclass$robpca_outlier
))

data_multiclass_clean <- multiclass_robpca %>%
  dplyr::filter(is.na(robpca_outlier) | robpca_outlier == FALSE) %>%
  dplyr::mutate(group_multiclass = factor(group_multiclass))

# Binary data are derived from the multiclass-cleaned dataset so both
# classification tasks use the same outlier-screened samples.
data_binary_clean <- data_multiclass_clean %>%
  dplyr::mutate(
    group_binary = factor(
      create_binary_labels(group_multiclass),
      levels = c("Others", "Multinational")
    )
  )

cat("\nMulticlass class distribution after outlier removal:\n")
print(table(data_multiclass_clean$group_multiclass))

cat("\nBinary class distribution after outlier removal:\n")
print(table(data_binary_clean$group_binary))


# ------------------------------------------------------------
# 8. LDA workflows
# ------------------------------------------------------------

multiclass_results <- run_lda_workflow(
  data = data_multiclass_clean,
  group_col = "group_multiclass",
  predictor_vars = inorganic_vars,
  model_label = "Inorganic multiclass LDA",
  output_prefix = "inorganic_multiclass_lda",
  seed = 123,
  training_proportion = TRAINING_PROPORTION
)

binary_results <- run_lda_workflow(
  data = data_binary_clean,
  group_col = "group_binary",
  predictor_vars = inorganic_vars,
  model_label = "Inorganic binary LDA",
  output_prefix = "inorganic_binary_lda",
  seed = 123,
  training_proportion = TRAINING_PROPORTION
)


# ------------------------------------------------------------
# 9. Robustness assessment
# ------------------------------------------------------------

if (RUN_ROBUSTNESS) {
  run_robustness_analysis(
    data = data_multiclass_clean,
    group_col = "group_multiclass",
    predictor_vars = inorganic_vars,
    output_prefix = "inorganic_multiclass_lda",
    n_repetitions = N_REPETITIONS,
    n_permutations = N_PERMUTATIONS,
    training_proportion = TRAINING_PROPORTION
  )

  run_robustness_analysis(
    data = data_binary_clean,
    group_col = "group_binary",
    predictor_vars = inorganic_vars,
    output_prefix = "inorganic_binary_lda",
    n_repetitions = N_REPETITIONS,
    n_permutations = N_PERMUTATIONS,
    training_proportion = TRAINING_PROPORTION
  )
}


# ------------------------------------------------------------
# 10. Session information
# ------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo_inorganic_model.txt")
)

cat("\nInorganic LDA workflow completed successfully.\n")
cat("Tables saved to: ", tables_dir, "\n")
cat("Figures saved to: ", figures_dir, "\n")
