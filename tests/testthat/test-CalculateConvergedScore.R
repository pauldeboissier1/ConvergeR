test_that("CalculateConvergedScore computes strictly numeric scores and correct columns", {
  # Charger le mini-objet PBMC truqué depuis le package
  seu_path <- system.file("extdata", "test_seu.rds", package = "ConvergeR")
  seu <- readRDS(seu_path)

  # Lancer ta fonction
  seu_res <- CalculateConvergedScore(
    seurat_obj = seu,
    principal_score = "scoreSeurat_Muscle",
    other_scores = "scoreSeurat_Fibrous",
    principal_name = "Muscle",
    other_name = "Fibrous"
  )

  # Vérifier que les colonnes ont bien été injectées dans les métadonnées
  expect_true("Converged_Score" %in% colnames(seu_res@meta.data))
  expect_true("Converged_Score_Direction" %in% colnames(seu_res@meta.data))

  # Vérifier que le calcul ne sort pas de valeurs aberrantes (ex: NA non justifiés)
  expect_type(seu_res@meta.data$Converged_Score, "double")
  expect_false(any(is.infinite(seu_res@meta.data$Converged_Score)))
})
