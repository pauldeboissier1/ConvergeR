library(Seurat)

# ==============================================================================
# SETUP: Creation of a mock Seurat mini-object for testing
# ==============================================================================
seu <- SeuratObject::pbmc_small
set.seed(42)

# Simulate module scores
seu$scoreSeurat_Muscle <- rnorm(ncol(seu), mean = 0.5, sd = 0.2)
seu$scoreSeurat_Fibrous <- rnorm(ncol(seu), mean = 0.3, sd = 0.2)

# Simulate clinical metadata
seu$age <- sample(20:70, ncol(seu), replace = TRUE)
seu$tumor_status <- sample(c("Primary", "Metastatic"), ncol(seu), replace = TRUE)
seu$treatment <- sample(c("Treated", "Untreated"), ncol(seu), replace = TRUE)
seu$SS <- paste0("SS", sample(1:5, ncol(seu), replace = TRUE))


# ==============================================================================
# TEST 1: Convergence score calculation
# ==============================================================================
test_that("CalculateConvergedScore generates correct metadata columns", {

  seu_res <- CalculateConvergedScore(
    seurat_obj = seu,
    principal_score = "scoreSeurat_Muscle",
    other_scores = "scoreSeurat_Fibrous",
    principal_name = "Muscle",
    other_name = "Fibrous",
    output_colname = "Converged_Score"
  )

  # Verify that the columns were successfully created
  expect_true("Converged_Score" %in% colnames(seu_res@meta.data))
  expect_true("Converged_Score_Direction" %in% colnames(seu_res@meta.data))
  expect_true("Converged_Score_Direction_Color" %in% colnames(seu_res@meta.data))

  # Verify that the score is strictly numeric without infinite values
  expect_type(seu_res@meta.data$Converged_Score, "double")
  expect_false(any(is.infinite(seu_res@meta.data$Converged_Score)))
})


# ==============================================================================
# TEST 2: Plotting functions (ggplot2)
# ==============================================================================
test_that("Plot functions return valid ggplot objects", {

  # Pre-calculate the score required for the plots
  seu_res <- CalculateConvergedScore(seu, "scoreSeurat_Muscle", "scoreSeurat_Fibrous")

  # Test the simple barplot
  p1 <- PlotConvergenceProportion(seu_res, x_var = "SS", fill_var = "Converged_Score_Direction")
  expect_s3_class(p1, "ggplot")

  # Test the crossed (faceted) barplot
  p2 <- PlotConvergenceCrossedProportion(seu_res, x_var = "SS", fill_var = "Converged_Score_Direction", facet_var = "tumor_status")
  expect_s3_class(p2, "ggplot")

  # Test the density plot
  p3 <- PlotConvergenceDensity(seu_res, x_var = "age", fill_var = "Converged_Score_Direction")
  expect_s3_class(p3, "ggplot")
})


# ==============================================================================
# TEST 3: Statistical functions
# ==============================================================================
test_that("Statistical tests return appropriate model objects", {

  seu_res <- CalculateConvergedScore(seu, "scoreSeurat_Muscle", "scoreSeurat_Fibrous")

  # Bivariate test (suppressMessages prevents console clutter during testing)
  res_biv <- suppressMessages(
    TestConvergenceScore(
      seurat_obj = seu_res,
      score_var = "Converged_Score",
      test_var = "age",
      level = "patient",
      patient_id_var = "SS"
    )
  )
  # htest is the standard R test class (wilcox, cor, etc.)
  expect_s3_class(res_biv, "htest")

  # Multivariate test
  res_multi <- suppressMessages(
    TestMultivariateConvergence(
      seurat_obj = seu_res,
      score_var = "Converged_Score",
      test_vars = c("age", "tumor_status"),
      level = "patient",
      patient_id_var = "SS"
    )
  )
  # Verify that the linear model ran successfully
  expect_s3_class(res_multi, "summary.lm")
})
