# ---------------------------------------------------------------------------
# Internal helpers: a single point of contact with the object class.
# Every exported function reads/writes cell metadata only through these,
# so Seurat and SingleCellExperiment (any SummarizedExperiment) are
# handled identically.
# ---------------------------------------------------------------------------

.is_sce <- function(x) methods::is(x, "SummarizedExperiment")
.is_seurat <- function(x) inherits(x, "Seurat")

.check_object <- function(x) {
  if (!(.is_sce(x) || .is_seurat(x))) {
    stop("`object` must be a SingleCellExperiment or a Seurat object.",
         call. = FALSE)
  }
  invisible(x)
}

.require_seurat <- function() {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("Package 'SeuratObject' is required to work with Seurat objects.",
         call. = FALSE)
  }
}

# Cell metadata as a data.frame (colData for SCE, meta.data for Seurat)
.get_meta <- function(x) {
  .check_object(x)
  if (.is_sce(x)) {
    as.data.frame(SummarizedExperiment::colData(x), optional = TRUE)
  } else {
    .require_seurat()
    x[[]]
  }
}

# Add (or overwrite) columns in the cell metadata.
# `new_cols` is a named list of vectors, each of length ncol(x).
.set_meta <- function(x, new_cols) {
  .check_object(x)
  if (.is_sce(x)) {
    cd <- SummarizedExperiment::colData(x)
    for (nm in names(new_cols)) {
      cd[[nm]] <- new_cols[[nm]]
    }
    SummarizedExperiment::colData(x) <- cd
    x
  } else {
    .require_seurat()
    md <- as.data.frame(new_cols, optional = TRUE, stringsAsFactors = FALSE)
    colnames(md) <- names(new_cols)
    rownames(md) <- colnames(x)
    SeuratObject::AddMetaData(x, metadata = md)
  }
}

.check_columns <- function(df, cols) {
  missing_cols <- setdiff(cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("The following columns are missing from the cell metadata: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

# Named color vector (direction label -> hex color) built from the metadata
.direction_colors <- function(df_meta, fill_var, color_var) {
  .check_columns(df_meta, c(fill_var, color_var))
  color_table <- unique(df_meta[, c(fill_var, color_var)])
  color_table <- color_table[stats::complete.cases(color_table), ]
  stats::setNames(color_table[[color_var]], color_table[[fill_var]])
}

# ---------------------------------------------------------------------------
# Exported functions
# ---------------------------------------------------------------------------

#' Calculate Converged Score
#'
#' Computes a "convergence score" comparing two groups of module
#' scores. Each group can be made of a single score or a combination of several
#' scores; the underlying methodology is identical either way.
#'
#' The function works on the cell-level metadata only (\code{colData} for a
#' \linkS4class{SingleCellExperiment}, \code{meta.data} for a Seurat object).
#' The object is never converted: it is returned with its original class and
#' all its content (assays, reduced dimensions, ...) untouched, with the new
#' columns added to the cell metadata.
#'
#' @param object A \linkS4class{SingleCellExperiment} (or any
#'   \linkS4class{SummarizedExperiment}) or a Seurat object. Its cell metadata
#'   must already contain the score columns.
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
#' @return The same object, with the convergence score columns added to its
#'   cell metadata.
#'
#' @references Whitfield HJ, Anderson ND, Burke C, Groot Koerkamp MJ, Parks C, Ogbonnah T, Wood Y, Piapi A, Robertson E, Watt E, White A, De Noon S, Kennedy J, Nagrecha R, Meister MT, Aladowicz E, Man YKS, Laspidea V, Trinh MK, Hodder A, Porter T, Lawrence JE, Tuck E, Nguyen T, Kelsey A, Flanagan AM, Hewitt R, Smeulders N, Slater O, Hutchinson JC, Sebire N, Shipley JM, Drost J, Straathof K, Behjati S. High-Risk Rhabdomyosarcomas Feature a Convergent Cell State. Cancer Res. 2026 Jul 14. doi: 10.1158/0008-5472.CAN-25-4403. PMID: 42446905.
#'
#' @examples
#' ## --- SingleCellExperiment ---
#' set.seed(42)
#' counts <- matrix(rpois(20000, 5), nrow = 200,
#'                  dimnames = list(paste0("gene", 1:200), paste0("cell", 1:100)))
#' sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
#'
#' # Simulated module scores (fast)
#' sce$Score_IFN <- rnorm(ncol(sce))
#' sce$Score_Inflam <- rnorm(ncol(sce))
#'
#' sce <- CalculateConvergedScore(
#'   object = sce,
#'   principal_score = "Score_IFN",
#'   other_scores = "Score_Inflam",
#'   principal_name = "Interferon",
#'   other_name = "Inflammatory",
#'   output_colname = "Converged_Immune",
#'   principal_color = "#9b59b6",
#'   other_color = "#e67e22"
#' )
#' head(SummarizedExperiment::colData(sce)[, c("Converged_Immune",
#'                                             "Converged_Immune_Direction")])
#'
#' ## --- Seurat (same function, same arguments) ---
#' if (requireNamespace("SeuratObject", quietly = TRUE)) {
#'   seu <- SeuratObject::CreateSeuratObject(counts = counts)
#'   seu$Score_IFN <- rnorm(ncol(seu))
#'   seu$Score_Inflam <- rnorm(ncol(seu))
#'   seu <- CalculateConvergedScore(
#'     object = seu,
#'     principal_score = "Score_IFN",
#'     other_scores = "Score_Inflam",
#'     output_colname = "Converged_Immune"
#'   )
#'   head(seu[[]][, c("Converged_Immune", "Converged_Immune_Direction")])
#' }
#'
#' \donttest{
#' # Full workflow with MSigDB signatures on a SingleCellExperiment
#' # (slower, requires scuttle and msigdbr)
#' if (requireNamespace("scuttle", quietly = TRUE) &&
#'     requireNamespace("msigdbr", quietly = TRUE)) {
#'   sce_path <- system.file("extdata", "pbmc3k_subset_sce.rds", package = "ConvergeR")
#'   sce <- readRDS(sce_path)
#'   sce <- scuttle::logNormCounts(sce)
#'   hs <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
#'   raw_group1 <- hs[hs$gs_name == "HALLMARK_INTERFERON_ALPHA_RESPONSE", ]$gene_symbol
#'   raw_group2 <- hs[hs$gs_name == "HALLMARK_INFLAMMATORY_RESPONSE", ]$gene_symbol
#'
#'   group1_genes <- intersect(raw_group1, rownames(sce))
#'   group2_genes <- intersect(raw_group2, rownames(sce))
#'
#'   # Mean log-normalized expression of each signature
#'   lc <- SingleCellExperiment::logcounts(sce)
#'   sce$Score_IFN <- colMeans(as.matrix(lc[group1_genes, , drop = FALSE]))
#'   sce$Score_Inflam <- colMeans(as.matrix(lc[group2_genes, , drop = FALSE]))
#'
#'   sce <- CalculateConvergedScore(
#'     object = sce,
#'     principal_score = "Score_IFN",
#'     other_scores = "Score_Inflam",
#'     output_colname = "Converged_Immune_MSigDB"
#'   )
#' }
#' }
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @import ggplot2
#' @import dplyr
#' @importFrom magrittr %>%
#' @import Seurat
#' @import SeuratObject
#' @importFrom stats as.formula complete.cases cor.test kruskal.test lm median na.omit pf sd setNames symnum wilcox.test
#' @export
CalculateConvergedScore <- function(object,
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

  meta <- .get_meta(object)

  # All new columns are collected here, then written once at the end.
  new_cols <- list()

  # Step 1: standardize (Z-score) every individual score, whether it
  # belongs to principal_score or other_scores.
  # unique() avoids recomputing a score twice if it happened to appear
  # in both groups.
  all_scores <- unique(c(principal_score, other_scores))
  .check_columns(meta, all_scores)

  for (sc in all_scores) {
    v <- meta[[sc]]
    if (!is.numeric(v)) {
      stop("Score '", sc, "' is not numeric.", call. = FALSE)
    }
    sd_sc <- stats::sd(v, na.rm = TRUE)
    if (is.na(sd_sc) || sd_sc == 0) {
      stop("Score '", sc, "' has zero (or undefined) variance and cannot be standardized.",
           call. = FALSE)
    }
    new_cols[[paste0(sc, "_Zscore")]] <- (v - mean(v, na.rm = TRUE)) / sd_sc
  }

  # Step 2: composite Z-score for each group (mean of Z-scores).
  # do.call(cbind, ...) always yields a matrix, even when a group has a
  # single score, so rowMeans() behaves the same for 1 or several scores.
  principal_zscore_cols <- paste0(principal_score, "_Zscore")
  other_zscore_cols     <- paste0(other_scores, "_Zscore")

  composite_principal <- rowMeans(
    do.call(cbind, new_cols[principal_zscore_cols]), na.rm = TRUE
  )
  composite_other <- rowMeans(
    do.call(cbind, new_cols[other_zscore_cols]), na.rm = TRUE
  )
  new_cols$Composite_Principal_ZScore <- composite_principal
  new_cols$Composite_Other_ZScore     <- composite_other

  # Step 3: rescale both composite Z-scores strictly to [0, 1]
  scaled_principal <- scales::rescale(composite_principal, to = c(0, 1))
  scaled_other     <- scales::rescale(composite_other, to = c(0, 1))
  new_cols$Scaled_Principal_ZScore <- scaled_principal
  new_cols$Scaled_Other_ZScore     <- scaled_other

  # Step 4: final convergence score (log2 ratio)
  # A constant of 0.001 is added to prevent log2(0)
  score <- log(scaled_principal + 0.001, base = 2) -
    log(scaled_other + 0.001, base = 2)
  new_cols[[output_colname]] <- score

  # Step 5: interpretation
  principal_label <- paste0("Convergence to ", principal_name, " (>0)")
  other_label     <- paste0("Convergence to ", other_name, " (<0)")
  direction <- ifelse(score > 0, principal_label, other_label)
  new_cols[[direction_colname]] <- direction

  # Step 6: matching hex color for each direction (one column per score,
  # named after direction_colname so parallel calls don't collide)
  new_cols[[color_colname]] <- ifelse(direction == principal_label,
                                      principal_color, other_color)

  .set_meta(object, new_cols)
}

#' Plot Convergence Proportion
#'
#' Draws a 100%-stacked barplot of cell proportions across categories of
#' x_var, colored by a convergence-direction column (fill_var). The color
#' mapping is extracted automatically from the cell metadata based on the
#' color column generated by \code{\link{CalculateConvergedScore}}.
#'
#' @param object A \linkS4class{SingleCellExperiment} or Seurat object whose
#'   cell metadata holds fill_var and its matching color column.
#' @param x_var Character. Name of the metadata column used on the x-axis (e.g. a patient identifier).
#' @param fill_var Character. Name of the metadata column used to fill the bars.
#' @param color_var Character or NULL. Name of the metadata column holding the hex color. Defaults to paste0(fill_var, "_Color").
#' @param title Character or NULL. Plot title. Defaults to x_var.
#'
#' @return A ggplot object.
#' @examples
#' set.seed(42)
#' counts <- matrix(rpois(20000, 5), nrow = 200,
#'                  dimnames = list(paste0("gene", 1:200), paste0("cell", 1:100)))
#' sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
#' sce$SS <- paste0("Patient_", sample(1:5, ncol(sce), replace = TRUE))
#' sce$Score_IFN <- rnorm(ncol(sce))
#' sce$Score_Inflam <- rnorm(ncol(sce))
#'
#' sce <- CalculateConvergedScore(
#'   object = sce,
#'   principal_score = "Score_IFN",
#'   other_scores = "Score_Inflam",
#'   output_colname = "Converged_Immune"
#' )
#'
#' p <- PlotConvergenceProportion(
#'   object = sce,
#'   x_var = "SS",
#'   fill_var = "Converged_Immune_Direction",
#'   title = "Immune Convergence by Patient"
#' )
#' print(p)
#' @importFrom ggplot2 ggplot aes geom_bar geom_hline scale_fill_manual
#'   theme_classic theme element_text labs
#' @importFrom rlang .data
#' @export
PlotConvergenceProportion <- function(object,
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

  df_meta <- .get_meta(object)
  .check_columns(df_meta, x_var)
  cols <- .direction_colors(df_meta, fill_var, color_var)

  ggplot(df_meta, aes(x = .data[[x_var]], fill = .data[[fill_var]])) +
    geom_bar(position = "fill") +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "black", linewidth = 1) +
    scale_fill_manual(values = cols) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title, x = NULL, y = "Cell Proportion")
}

#' Plot Convergence Crossed Proportion
#'
#' Draws a 100%-stacked barplot of cell proportions across categories of
#' x_var, colored by a convergence-direction column, and separated (faceted)
#' by a secondary categorical variable (facet_var).
#'
#' @param object A \linkS4class{SingleCellExperiment} or Seurat object whose
#'   cell metadata holds x_var, fill_var, facet_var, and the color column.
#' @param x_var Character. Name of the metadata column used on the x-axis.
#' @param fill_var Character. Name of the metadata column used to fill the bars.
#' @param facet_var Character. Name of the metadata column used to facet the plot (e.g., "tumor_status").
#' @param color_var Character or NULL. Name of the metadata column holding the hex color.
#' @param title Character or NULL. Plot title.
#' @param x_label Character or NULL. X-axis label.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' set.seed(42)
#' counts <- matrix(rpois(20000, 5), nrow = 200,
#'                  dimnames = list(paste0("gene", 1:200), paste0("cell", 1:100)))
#' sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
#' sce$SS <- paste0("Patient_", sample(1:5, ncol(sce), replace = TRUE))
#' sce$tumor_status <- sample(c("Primary", "Metastatic"), ncol(sce), replace = TRUE)
#' sce$Score_IFN <- rnorm(ncol(sce))
#' sce$Score_Inflam <- rnorm(ncol(sce))
#'
#' sce <- CalculateConvergedScore(
#'   object = sce,
#'   principal_score = "Score_IFN",
#'   other_scores = "Score_Inflam",
#'   output_colname = "Converged_Immune"
#' )
#'
#' b_cross <- PlotConvergenceCrossedProportion(
#'   object = sce,
#'   x_var = "SS",
#'   fill_var = "Converged_Immune_Direction",
#'   facet_var = "tumor_status",
#'   title = "Immune Convergence split by Tumor Status",
#'   x_label = "Patient ID"
#' )
#' print(b_cross)
#' @importFrom ggplot2 ggplot aes geom_bar geom_hline scale_fill_manual
#'   facet_grid vars theme_classic theme element_text labs
#' @importFrom rlang .data
PlotConvergenceCrossedProportion <- function(object,
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

  df_meta <- .get_meta(object)
  .check_columns(df_meta, c(x_var, facet_var))
  cols <- .direction_colors(df_meta, fill_var, color_var)

  ggplot(df_meta, aes(x = .data[[x_var]], fill = .data[[fill_var]])) +
    geom_bar(position = "fill") +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "black", linewidth = 1) +
    scale_fill_manual(values = cols) +
    facet_grid(cols = vars(.data[[facet_var]]), scales = "free_x", space = "free_x") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title, x = x_label, y = "Cell Proportion")
}

#' Plot Convergence Density
#'
#' Draws overlapping density curves of a continuous variable (x_var),
#' split by a convergence-direction column (fill_var).
#'
#' @param object A \linkS4class{SingleCellExperiment} or Seurat object whose
#'   cell metadata holds x_var, fill_var, and its matching color column.
#' @param x_var Character. Name of the metadata column holding the continuous variable (e.g. "age").
#' @param fill_var Character. Name of the metadata column used to fill the density curves.
#' @param color_var Character or NULL. Name of the metadata column holding the hex color.
#' @param title Character or NULL. Plot title.
#' @param x_label Character or NULL. X-axis label. Defaults to x_var.
#' @param alpha Numeric. Transparency of the density fill, between 0 and 1 (default: 0.6).
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' set.seed(42)
#' counts <- matrix(rpois(20000, 5), nrow = 200,
#'                  dimnames = list(paste0("gene", 1:200), paste0("cell", 1:100)))
#' sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
#' sce$age <- sample(20:70, ncol(sce), replace = TRUE)
#' sce$Score_IFN <- rnorm(ncol(sce))
#' sce$Score_Inflam <- rnorm(ncol(sce))
#'
#' sce <- CalculateConvergedScore(
#'   object = sce,
#'   principal_score = "Score_IFN",
#'   other_scores = "Score_Inflam",
#'   output_colname = "Converged_Immune"
#' )
#'
#' PlotConvergenceDensity(
#'   object = sce,
#'   x_var = "age",
#'   fill_var = "Converged_Immune_Direction",
#'   x_label = "Patient Age"
#' )
#' @importFrom ggplot2 ggplot aes geom_density scale_fill_manual
#'   theme_classic theme labs
#' @importFrom rlang .data
PlotConvergenceDensity <- function(object,
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

  df_meta <- .get_meta(object)
  .check_columns(df_meta, x_var)
  cols <- .direction_colors(df_meta, fill_var, color_var)

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
#' @param object A \linkS4class{SingleCellExperiment} or Seurat object whose
#'   cell metadata holds score_var and test_var.
#' @param score_var Character. Name of the metadata column holding the numeric score.
#' @param test_var Character. Name of the metadata column to test the score against.
#' @param level Character. Either "cell" (default) or "patient" (pseudobulk aggregation).
#' @param patient_id_var Character. Name of the metadata column identifying patients (required if level = "patient").
#'
#' @return The underlying htest object (invisibly). A message is printed summarizing the results.
#' @importFrom stats median na.omit complete.cases cor.test wilcox.test kruskal.test
#' @importFrom dplyr group_by summarise transmute
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' set.seed(42)
#' counts <- matrix(rpois(20000, 5), nrow = 200,
#'                  dimnames = list(paste0("gene", 1:200), paste0("cell", 1:100)))
#' sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
#' pat <- paste0("Patient_", 1:5)
#' sce$SS <- sample(pat, ncol(sce), replace = TRUE)
#' # age is a patient-level covariate (constant within a patient)
#' sce$age <- c(25, 33, 47, 52, 68)[match(sce$SS, pat)]
#' sce$Score_IFN <- rnorm(ncol(sce))
#' sce$Score_Inflam <- rnorm(ncol(sce))
#'
#' sce <- CalculateConvergedScore(
#'   object = sce,
#'   principal_score = "Score_IFN",
#'   other_scores = "Score_Inflam",
#'   output_colname = "Converged_Immune"
#' )
#'
#' TestConvergenceScore(
#'   object = sce,
#'   score_var = "Converged_Immune",
#'   test_var = "age",
#'   level = "patient",
#'   patient_id_var = "SS"
#' )
TestConvergenceScore <- function(object,
                                 score_var,
                                 test_var,
                                 level = c("cell", "patient"),
                                 patient_id_var = NULL) {

  level <- match.arg(level)

  if (level == "patient" && is.null(patient_id_var)) {
    stop("patient_id_var must be provided when level = \"patient\" (needed to aggregate one row per patient).")
  }

  df_meta <- .get_meta(object)
  .check_columns(df_meta, c(score_var, test_var,
                            if (level == "patient") patient_id_var))

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
#' @param object A \linkS4class{SingleCellExperiment} or Seurat object whose
#'   cell metadata holds the variables.
#' @param score_var Character. Name of the metadata column holding the numeric score.
#' @param test_vars Character vector. Names of the metadata columns to include in the model.
#' @param level Character. Either "cell" (default) or "patient" (pseudobulk aggregation).
#' @param patient_id_var Character. Name of the metadata column identifying patients (required if level = "patient").
#'
#' @return The summary.lm() object (invisibly). A detailed message is printed summarizing the model.
#' @examples
#' set.seed(42)
#' counts <- matrix(rpois(20000, 5), nrow = 200,
#'                  dimnames = list(paste0("gene", 1:200), paste0("cell", 1:100)))
#' sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
#' pat <- paste0("Patient_", 1:8)
#' sce$SS <- sample(pat, ncol(sce), replace = TRUE)
#' # patient-level covariates (constant within a patient)
#' sce$age <- c(25, 33, 47, 52, 68, 41, 59, 36)[match(sce$SS, pat)]
#' sce$tumor_status <- rep(c("Primary", "Metastatic"), 4)[match(sce$SS, pat)]
#' sce$Score_IFN <- rnorm(ncol(sce))
#' sce$Score_Inflam <- rnorm(ncol(sce))
#'
#' sce <- CalculateConvergedScore(
#'   object = sce,
#'   principal_score = "Score_IFN",
#'   other_scores = "Score_Inflam",
#'   output_colname = "Converged_Immune"
#' )
#'
#' TestMultivariateConvergence(
#'   object = sce,
#'   score_var = "Converged_Immune",
#'   test_vars = c("age", "tumor_status"),
#'   level = "patient",
#'   patient_id_var = "SS"
#' )
#' @importFrom stats median na.omit complete.cases lm as.formula pf symnum
#' @importFrom dplyr group_by summarise across all_of
#' @importFrom rlang .data
#' @export
TestMultivariateConvergence <- function(object,
                                        score_var,
                                        test_vars,
                                        level = c("cell", "patient"),
                                        patient_id_var = NULL) {

  level <- match.arg(level)

  if (level == "patient" && is.null(patient_id_var)) {
    stop("patient_id_var must be provided when level = \"patient\" (needed to aggregate one row per patient).")
  }

  df_meta <- .get_meta(object)
  .check_columns(df_meta, c(score_var, test_vars,
                            if (level == "patient") patient_id_var))

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
    df_test <- df_meta[, c(score_var, test_vars), drop = FALSE]
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
