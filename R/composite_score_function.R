#' Calculate Converged Score
#'
#' Computes a "convergence score" comparing two groups of Seurat module
#' scores. Each group can be made of a single score or a combination of several
#' scores; the underlying methodology is identical either way.
#'
#' @param seurat_obj Seurat object. Its meta.data must already contain the score columns.
#' @param principal_score Character. Name of one score, OR a character vector naming several scores.
#' @param other_scores Character. Same as principal_score, but for the group compared against.
#' @param principal_name Character. Human-readable label (default: "Principal").
#' @param other_name Character. Human-readable label (default: "Other").
#' @param output_colname Character. Name of the column for the score (default: "Converged_Score").
#' @param direction_colname Character or NULL. Name of the interpretation column.
#' @param color_colname Character or NULL. Name of the hex color column.
#' @param principal_color Character. Hex color for principal group (default: "#2ecc71").
#' @param other_color Character. Hex color for other group (default: "#e74c3c").
#'
#' @return The same Seurat object, with the convergence score columns added to meta.data.
#'
#' @references Whitfield HJ, Anderson ND, Burke C, Groot Koerkamp MJ, Parks C, Ogbonnah T, Wood Y, Piapi A, Robertson E, Watt E, White A, De Noon S, Kennedy J, Nagrecha R, Meister MT, Aladowicz E, Man YKS, Laspidea V, Trinh MK, Hodder A, Porter T, Lawrence JE, Tuck E, Nguyen T, Kelsey A, Flanagan AM, Hewitt R, Smeulders N, Slater O, Hutchinson JC, Sebire N, Shipley JM, Drost J, Straathof K, Behjati S. High-Risk Rhabdomyosarcomas Feature a Convergent Cell State. Cancer Res. 2026 Jul 14. doi: 10.1158/0008-5472.CAN-25-4403. PMID: 42446905.
#'
#' @examples
#' # 1. Environment Setup and Data Preparation
#' seu_path <- system.file("extdata", "pbmc3k_subset.rds", package = "ConvergeR")
#' seu <- readRDS(seu_path)
#' seu <- Seurat::NormalizeData(seu, verbose = FALSE)
#'
#' set.seed(42)
#' seu$age <- sample(20:70, ncol(seu), replace = TRUE)
#' seu$tumor_status <- sample(c("Primary", "Metastatic"), ncol(seu), replace = TRUE)
#' seu$treatment <- sample(c("Treated", "Untreated"), ncol(seu), replace = TRUE)
#' seu$SS <- paste0("Patient_", sample(1:5, ncol(seu), replace = TRUE))
#'
#' # 2. Importing MSigDB Signatures and Scoring
#' if (requireNamespace("msigdbr", quietly = TRUE)) {
#'   hs <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
#'   raw_group1 <- hs[hs$gs_name == "HALLMARK_INTERFERON_ALPHA_RESPONSE", ]$gene_symbol
#'   raw_group2 <- hs[hs$gs_name == "HALLMARK_INFLAMMATORY_RESPONSE", ]$gene_symbol
#'
#'   group1_genes <- intersect(raw_group1, rownames(seu))
#'   group2_genes <- intersect(raw_group2, rownames(seu))
#'
#'   seu <- Seurat::AddModuleScore(seu, features = list(group1_genes), name = "Score_IFN_")
#'   seu <- Seurat::AddModuleScore(seu, features = list(group2_genes), name = "Score_Inflam_")
#'   seu$Score_IFN <- seu$Score_IFN_1
#'   seu$Score_Inflam <- seu$Score_Inflam_1
#'
#'   # 3. Computing the Convergence Score
#'   seu <- CalculateConvergedScore(
#'     seurat_obj = seu,
#'     principal_score = "Score_IFN",
#'     other_scores = "Score_Inflam",
#'     principal_name = "Interferon",
#'     other_name = "Inflammatory",
#'     output_colname = "Converged_Immune",
#'     principal_color = "#9b59b6",
#'     other_color = "#e67e22"
#'   )
#' }
#' @import ggplot2
#' @import dplyr
#' @importFrom magrittr %>%
#' @import Seurat
#' @import SeuratObject
#' @importFrom stats as.formula complete.cases cor.test kruskal.test lm median na.omit pf sd setNames symnum wilcox.test
#' @export
CalculateConvergedScore <- function(seurat_obj,
                                    principal_score,
                                    other_scores,
                                    principal_name = "Principal",
                                    other_name = "Other",
                                    output_colname = "Converged_Score",
                                    direction_colname = NULL,
                                    color_colname = NULL,
                                    principal_color = "#2ecc71",
                                    other_color = "#e74c3c") {

  if (is.null(direction_colname)) {
    direction_colname <- paste0(output_colname, "_Direction")
  }
  if (is.null(color_colname)) {
    color_colname <- paste0(direction_colname, "_Color")
  }

  # Step 1: standardize (Z-score) every individual score, whether it
  # belongs to principal_score or other_scores.
  # unique() avoids recomputing a score twice if it happened to appear
  # in both groups.
  all_scores <- unique(c(principal_score, other_scores))

  for (sc in all_scores) {
    mean_sc <- mean(seurat_obj@meta.data[[sc]], na.rm = TRUE)
    sd_sc   <- sd(seurat_obj@meta.data[[sc]], na.rm = TRUE)

    zscore_colname <- paste0(sc, "_Zscore")
    seurat_obj@meta.data[[zscore_colname]] <- (seurat_obj@meta.data[[sc]] - mean_sc) / sd_sc
  }

  # Step 2: composite Z-score for each group (mean of Z-scores).
  # drop = FALSE is essential here: it keeps the selection as a
  # data.frame even when a group has only one score, so rowMeans()
  # works the same way whether the group has 1 or several scores.
  principal_zscore_cols <- paste0(principal_score, "_Zscore")
  other_zscore_cols     <- paste0(other_scores, "_Zscore")

  seurat_obj@meta.data$Composite_Principal_ZScore <- rowMeans(
    seurat_obj@meta.data[, principal_zscore_cols, drop = FALSE], na.rm = TRUE
  )
  seurat_obj@meta.data$Composite_Other_ZScore <- rowMeans(
    seurat_obj@meta.data[, other_zscore_cols, drop = FALSE], na.rm = TRUE
  )

  # Step 3: rescale both composite Z-scores strictly to [0, 1]
  seurat_obj@meta.data$Scaled_Principal_ZScore <- scales::rescale(
    seurat_obj@meta.data$Composite_Principal_ZScore, to = c(0, 1)
  )
  seurat_obj@meta.data$Scaled_Other_ZScore <- scales::rescale(
    seurat_obj@meta.data$Composite_Other_ZScore, to = c(0, 1)
  )

  # Step 4: final convergence score (log2 ratio)
  # A constant of 0.001 is added to prevent log2(0)
  seurat_obj@meta.data[[output_colname]] <-
    log(seurat_obj@meta.data$Scaled_Principal_ZScore + 0.001, base = 2) -
    log(seurat_obj@meta.data$Scaled_Other_ZScore + 0.001, base = 2)

  # Step 5: interpretation, computed directly inside the function
  seurat_obj@meta.data[[direction_colname]] <- ifelse(
    seurat_obj@meta.data[[output_colname]] > 0,
    paste0("Convergence to ", principal_name, " (>0)"),
    paste0("Convergence to ", other_name, " (<0)")
  )

  # Step 6: matching hex color for each direction (one column per score,
  # named after direction_colname so parallel calls don't collide)
  seurat_obj@meta.data[[color_colname]] <- ifelse(
    seurat_obj@meta.data[[direction_colname]] == paste0("Convergence to ", principal_name, " (>0)"),
    principal_color,
    other_color
  )

  return(seurat_obj)
}

#' Plot Convergence Proportion
#'
#' Draws a 100%-stacked barplot of cell proportions across categories of
#' x_var, colored by a convergence-direction column (fill_var). The color
#' mapping is extracted automatically from meta.data based on the color column
#' generated by CalculateConvergedScore().
#'
#' @param seurat_obj Seurat object whose meta.data holds both fill_var and its matching color column.
#' @param x_var Character. Name of the meta.data column used on the x-axis (e.g. a patient identifier).
#' @param fill_var Character. Name of the meta.data column used to fill the bars.
#' @param color_var Character or NULL. Name of the meta.data column holding the hex color. Defaults to paste0(fill_var, "_Color").
#' @param title Character or NULL. Plot title. Defaults to x_var.
#'
#' @return A ggplot object.
#' @examples
#' seu_path <- system.file("extdata", "pbmc3k_subset.rds", package = "ConvergeR")
#' seu <- readRDS(seu_path)
#' set.seed(42)
#' seu$SS <- paste0("Patient_", sample(1:5, ncol(seu), replace = TRUE))
#'
#' # Simulate fast scores to ensure example runs efficiently
#' seu$Score_IFN <- rnorm(ncol(seu))
#' seu$Score_Inflam <- rnorm(ncol(seu))
#'
#' seu <- CalculateConvergedScore(
#'   seurat_obj = seu,
#'   principal_score = "Score_IFN",
#'   other_scores = "Score_Inflam",
#'   output_colname = "Converged_Immune"
#' )
#'
#' p <- PlotConvergenceProportion(
#'   seurat_obj = seu,
#'   x_var = "SS",
#'   fill_var = "Converged_Immune_Direction",
#'   title = "Immune Convergence by Patient"
#' )
#' print(p)
#' @export
PlotConvergenceProportion <- function(seurat_obj,
                                      x_var,
                                      fill_var,
                                      color_var = NULL,
                                      title = NULL) {

  if (is.null(color_var)) {
    color_var <- paste0(fill_var, "_Color")
  }
  if (is.null(title)) {
    title <- x_var
  }

  df_meta <- seurat_obj@meta.data

  # Automatically build the named color vector from meta.data: one
  # unique (direction label -> hex color) pair per category.
  color_table <- unique(df_meta[, c(fill_var, color_var)])
  cols <- setNames(color_table[[color_var]], color_table[[fill_var]])

  ggplot(df_meta, aes(x = .data[[x_var]], fill = .data[[fill_var]])) +
    geom_bar(position = "fill") +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "black", linewidth = 1) +
    scale_fill_manual(values = cols) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title, x = NULL, y = "Proportion de cellules")
}

#' Plot Convergence Crossed Proportion
#'
#' Draws a 100%-stacked barplot of cell proportions across categories of
#' x_var, colored by a convergence-direction column, and separated (faceted)
#' by a secondary categorical variable (facet_var).
#'
#' @param seurat_obj Seurat object whose meta.data holds x_var, fill_var, facet_var, and the color column.
#' @param x_var Character. Name of the meta.data column used on the x-axis.
#' @param fill_var Character. Name of the meta.data column used to fill the bars.
#' @param facet_var Character. Name of the meta.data column used to facet the plot (e.g., "tumor_status").
#' @param color_var Character or NULL. Name of the meta.data column holding the hex color.
#' @param title Character or NULL. Plot title.
#' @param x_label Character or NULL. X-axis label.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' seu_path <- system.file("extdata", "pbmc3k_subset.rds", package = "ConvergeR")
#' seu <- readRDS(seu_path)
#' set.seed(42)
#' seu$SS <- paste0("Patient_", sample(1:5, ncol(seu), replace = TRUE))
#' seu$tumor_status <- sample(c("Primary", "Metastatic"), ncol(seu), replace = TRUE)
#' seu$Score_IFN <- rnorm(ncol(seu))
#' seu$Score_Inflam <- rnorm(ncol(seu))
#'
#' seu <- CalculateConvergedScore(
#'   seurat_obj = seu,
#'   principal_score = "Score_IFN",
#'   other_scores = "Score_Inflam",
#'   output_colname = "Converged_Immune"
#' )
#'
#' b_cross <- PlotConvergenceCrossedProportion(
#'   seurat_obj = seu,
#'   x_var = "SS",
#'   fill_var = "Converged_Immune_Direction",
#'   facet_var = "tumor_status",
#'   title = "Immune Convergence split by Tumor Status",
#'   x_label = "Patient ID"
#' )
#' print(b_cross)
PlotConvergenceCrossedProportion <- function(seurat_obj,
                                             x_var,
                                             fill_var,
                                             facet_var,
                                             color_var = NULL,
                                             title = NULL,
                                             x_label = NULL) {

  if (is.null(color_var)) {
    color_var <- paste0(fill_var, "_Color")
  }
  if (is.null(title)) {
    title <- paste0("Cell Proportion across ", x_var, " split by ", facet_var)
  }
  if (is.null(x_label)) {
    x_label <- x_var
  }

  df_meta <- seurat_obj@meta.data

  # Automatically build the named color vector from meta.data: one
  # unique (direction label -> hex color) pair per category.
  color_table <- unique(df_meta[, c(fill_var, color_var)])
  cols <- setNames(color_table[[color_var]], color_table[[fill_var]])

  ggplot(df_meta, aes(x = .data[[x_var]], fill = .data[[fill_var]])) +
    geom_bar(position = "fill") +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "black", linewidth = 1) +
    scale_fill_manual(values = cols) +
    facet_grid(as.formula(paste("~", facet_var)), scales = "free_x", space = "free_x") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title, x = x_label, y = "Cell Proportion")
}

#' Plot Convergence Density
#'
#' Draws overlapping density curves of a continuous variable (x_var),
#' split by a convergence-direction column (fill_var).
#'
#' @param seurat_obj Seurat object whose meta.data holds x_var, fill_var, and its matching color column.
#' @param x_var Character. Name of the meta.data column holding the continuous variable (e.g. "age").
#' @param fill_var Character. Name of the meta.data column used to fill the density curves.
#' @param color_var Character or NULL. Name of the meta.data column holding the hex color.
#' @param title Character or NULL. Plot title.
#' @param x_label Character or NULL. X-axis label. Defaults to x_var.
#' @param alpha Numeric. Transparency of the density fill, between 0 and 1 (default: 0.6).
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' seu_path <- system.file("extdata", "pbmc3k_subset.rds", package = "ConvergeR")
#' seu <- readRDS(seu_path)
#' set.seed(42)
#' seu$age <- sample(20:70, ncol(seu), replace = TRUE)
#' seu$Score_IFN <- rnorm(ncol(seu))
#' seu$Score_Inflam <- rnorm(ncol(seu))
#'
#' seu <- CalculateConvergedScore(
#'   seurat_obj = seu,
#'   principal_score = "Score_IFN",
#'   other_scores = "Score_Inflam",
#'   output_colname = "Converged_Immune"
#' )
#'
#' PlotConvergenceDensity(
#'   seurat_obj = seu,
#'   x_var = "age",
#'   fill_var = "Converged_Immune_Direction",
#'   x_label = "Patient Age"
#' )
PlotConvergenceDensity <- function(seurat_obj,
                                   x_var,
                                   fill_var,
                                   color_var = NULL,
                                   title = NULL,
                                   x_label = NULL,
                                   alpha = 0.6) {

  if (is.null(color_var)) {
    color_var <- paste0(fill_var, "_Color")
  }
  if (is.null(title)) {
    title <- paste0(x_var, " Distribution by Convergence Direction")
  }
  if (is.null(x_label)) {
    x_label <- x_var
  }

  df_meta <- seurat_obj@meta.data

  # Automatically build the named color vector from meta.data: one
  # unique (direction label -> hex color) pair per category.
  color_table <- unique(df_meta[, c(fill_var, color_var)])
  cols <- setNames(color_table[[color_var]], color_table[[fill_var]])

  ggplot(df_meta, aes(x = .data[[x_var]], fill = .data[[fill_var]])) +
    geom_density(alpha = alpha) +
    scale_fill_manual(values = cols) +
    theme_classic() +
    labs(title = title, x = x_label, y = "Cell Density") +
    theme(legend.position = "bottom")
}

#' Test Convergence Score
#'
#' Runs a statistical test between a numeric convergence score and another variable,
#' automatically choosing the appropriate test (Spearman, Wilcoxon, or Kruskal-Wallis).
#'
#' @param seurat_obj Seurat object whose meta.data holds score_var and test_var.
#' @param score_var Character. Name of the meta.data column holding the numeric score.
#' @param test_var Character. Name of the meta.data column to test the score against.
#' @param level Character. Either "cell" (default) or "patient" (pseudobulk aggregation).
#' @param patient_id_var Character. Name of the meta.data column identifying patients (required if level = "patient").
#'
#' @return The underlying htest object (invisibly). A message is printed summarizing the results.
#' @importFrom stats sd median na.omit complete.cases cor.test wilcox.test kruskal.test lm as.formula pf symnum setNames
#' @export
#'
#' @examples
#' seu_path <- system.file("extdata", "pbmc3k_subset.rds", package = "ConvergeR")
#' seu <- readRDS(seu_path)
#' set.seed(42)
#' seu$SS <- paste0("Patient_", sample(1:5, ncol(seu), replace = TRUE))
#' seu$age <- sample(20:70, ncol(seu), replace = TRUE)
#' seu$Score_IFN <- rnorm(ncol(seu))
#' seu$Score_Inflam <- rnorm(ncol(seu))
#'
#' seu <- CalculateConvergedScore(
#'   seurat_obj = seu,
#'   principal_score = "Score_IFN",
#'   other_scores = "Score_Inflam",
#'   output_colname = "Converged_Immune"
#' )
#'
#' TestConvergenceScore(
#'   seurat_obj = seu,
#'   score_var = "Converged_Immune",
#'   test_var = "age",
#'   level = "patient",
#'   patient_id_var = "SS"
#' )
TestConvergenceScore <- function(seurat_obj,
                                 score_var,
                                 test_var,
                                 level = c("cell", "patient"),
                                 patient_id_var = NULL) {

  level <- match.arg(level)

  if (level == "patient" && is.null(patient_id_var)) {
    stop("patient_id_var must be provided when level = \"patient\" (needed to aggregate one row per patient).")
  }

  df_meta <- seurat_obj@meta.data

  if (level == "patient") {
    # Pseudobulk aggregation: one row per patient, correcting for
    # pseudoreplication. The convergence score is summarized by its
    # median per patient; test_var is assumed to be a patient-level
    # covariate (constant across a patient's cells), so its first
    # non-missing value is kept.
    df_test <- df_meta %>%
      group_by(.data[[patient_id_var]]) %>%
      summarise(
        .score = median(.data[[score_var]], na.rm = TRUE),
        .var = dplyr::first(na.omit(.data[[test_var]])),
        .groups = "drop"
      )
  } else {
    df_test <- df_meta %>%
      transmute(.score = .data[[score_var]], .var = .data[[test_var]])
  }

  df_test <- df_test[complete.cases(df_test[, c(".score", ".var")]), ]

  # Automatically pick the appropriate test based on the type of test_var
  if (is.numeric(df_test$.var)) {
    test_type <- "Spearman correlation"
    result <- cor.test(df_test$.score, df_test$.var, method = "spearman")
  } else {
    grp <- droplevels(as.factor(df_test$.var))
    n_levels <- nlevels(grp)

    if (n_levels < 2) {
      stop("test_var has fewer than 2 groups after removing missing values; a comparison test cannot be run.")
    } else if (n_levels == 2) {
      test_type <- "Wilcoxon rank-sum test"
      result <- wilcox.test(.score ~ .var, data = df_test)
    } else {
      test_type <- "Kruskal-Wallis test"
      result <- kruskal.test(.score ~ .var, data = df_test)
    }
  }

  message(
    test_type,
    " | level: ", ifelse(level == "patient", "patient (pseudobulk)", "cell"),
    " | N = ", nrow(df_test),
    " | testing '", score_var, "' by '", test_var, "'",
    " | p = ", format.pval(result$p.value, digits = 3)
  )

  return(invisible(result))
}

#' Test Multivariate Convergence
#'
#' Fits a multivariate linear model (lm) to test a numeric convergence score
#' against multiple covariates simultaneously, correcting for confounding factors.
#'
#' @param seurat_obj Seurat object whose meta.data holds the variables.
#' @param score_var Character. Name of the meta.data column holding the numeric score.
#' @param test_vars Character vector. Names of the meta.data columns to include in the model.
#' @param level Character. Either "cell" (default) or "patient" (pseudobulk aggregation).
#' @param patient_id_var Character. Name of the meta.data column identifying patients (required if level = "patient").
#'
#' @return The summary.lm() object (invisibly). A detailed message is printed summarizing the model.
#' @examples
#' seu_path <- system.file("extdata", "pbmc3k_subset.rds", package = "ConvergeR")
#' seu <- readRDS(seu_path)
#' set.seed(42)
#' seu$SS <- paste0("Patient_", sample(1:5, ncol(seu), replace = TRUE))
#' seu$age <- sample(20:70, ncol(seu), replace = TRUE)
#' seu$tumor_status <- sample(c("Primary", "Metastatic"), ncol(seu), replace = TRUE)
#' seu$treatment <- sample(c("Treated", "Untreated"), ncol(seu), replace = TRUE)
#' seu$Score_IFN <- rnorm(ncol(seu))
#' seu$Score_Inflam <- rnorm(ncol(seu))
#'
#' seu <- CalculateConvergedScore(
#'   seurat_obj = seu,
#'   principal_score = "Score_IFN",
#'   other_scores = "Score_Inflam",
#'   output_colname = "Converged_Immune"
#' )
#'
#' TestMultivariateConvergence(
#'   seurat_obj = seu,
#'   score_var = "Converged_Immune",
#'   test_vars = c("age", "tumor_status", "treatment"),
#'   level = "patient",
#'   patient_id_var = "SS"
#' )
#' @importFrom stats sd median na.omit complete.cases cor.test wilcox.test kruskal.test lm as.formula pf symnum setNames
#' @export
TestMultivariateConvergence <- function(seurat_obj,
                                        score_var,
                                        test_vars,
                                        level = c("cell", "patient"),
                                        patient_id_var = NULL) {

  level <- match.arg(level)

  if (level == "patient" && is.null(patient_id_var)) {
    stop("patient_id_var must be provided when level = \"patient\" (needed to aggregate one row per patient).")
  }

  df_meta <- seurat_obj@meta.data

  # Check if all columns exist
  missing_cols <- setdiff(c(score_var, test_vars, if(level=="patient") patient_id_var), colnames(df_meta))
  if (length(missing_cols) > 0) {
    stop("The following columns are missing from meta.data: ", paste(missing_cols, collapse = ", "))
  }

  if (level == "patient") {
    # Pseudobulk aggregation
    df_test <- df_meta %>%
      group_by(.data[[patient_id_var]]) %>%
      summarise(
        .score = median(.data[[score_var]], na.rm = TRUE),
        # Extract the first non-missing value for each covariate
        across(all_of(test_vars), ~ dplyr::first(na.omit(.x))),
        .groups = "drop"
      )
  } else {
    # Cell-level formatting
    df_test <- df_meta[, c(score_var, test_vars)]
    colnames(df_test)[1] <- ".score"
  }

  # Remove rows with NAs in any of the tested variables to ensure a clean model
  df_test <- df_test[complete.cases(df_test), ]

  # Build the formula dynamically: .score ~ var1 + var2 + ...
  form_str <- paste(".score ~", paste(test_vars, collapse = " + "))

  # Run the linear model
  mod <- lm(as.formula(form_str), data = df_test)
  mod_summary <- summary(mod)

  # Extract global p-value of the model (F-statistic)
  f_stat <- mod_summary$fstatistic
  if (!is.null(f_stat)) {
    global_p <- pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)
  } else {
    global_p <- NA
  }

  # Format the output message
  message("======================================================")
  message("Multivariate Linear Model | level: ", ifelse(level == "patient", "patient (pseudobulk)", "cell"), " | N = ", nrow(df_test))
  message("Formula: ", score_var, " ~ ", paste(test_vars, collapse = " + "))
  message("Global Model R-squared: ", round(mod_summary$r.squared, 3),
          " | Global p-value: ", format.pval(global_p, digits = 3))
  message("------------------------------------------------------")
  message("Independent Variable Effects:")

  # Loop through coefficients to print interpretation
  coef_table <- mod_summary$coefficients

  for (row_name in rownames(coef_table)) {
    if (row_name == "(Intercept)") next

    estimate <- coef_table[row_name, "Estimate"]
    p_val <- coef_table[row_name, "Pr(>|t|)"]

    # Assign significance stars for quick reading
    stars <- symnum(p_val, corr = FALSE, na = FALSE,
                    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                    symbols = c("***", "**", "*", ".", " "))

    message(sprintf("  - %-25s : p = %-8s %s (Estimate = %6.3f)",
                    row_name,
                    format.pval(p_val, digits = 3),
                    trimws(stars),
                    estimate))
  }
  message("======================================================")

  return(invisible(mod_summary))
}
