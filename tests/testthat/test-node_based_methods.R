data(toy_pathway, toy_gene_DEdata, toy_mir_DEdata, toy_gene_mat, toy_mir_mat,
     toy_background_MGI, toy_background_GGI, package = "miREA")

test_that("TG_ORA and TG_Score return enrichment results for the signal pathway", {
  input_data <- get_all_input_data(
    methods = c("TG_ORA", "TG_Score"), pathway = toy_pathway,
    mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
    background_MGI = toy_background_MGI
  )

  ora <- TG_ORA(TG = input_data$data$TG_ORA, background = unique(toy_pathway$gene),
                pathway = input_data$pathway$TG, pvalueCutoff = 1)
  expect_s3_class(ora, "data.frame")
  expect_true("REPLICATIVE IMMORTALITY" %in% ora$pathway)

  score <- TG_Score(geneList = input_data$data$TG_Score, pathway = input_data$pathway$TG,
                     pvalueType = "std", pvalueCutoff = 1)
  expect_s3_class(score, "data.frame")
  expect_true("REPLICATIVE IMMORTALITY" %in% score$pathway)
})

test_that("TG_ORA rejects an empty TG vector", {
  expect_error(TG_ORA(TG = character(0), pathway = toy_pathway), "TG is empty")
})

test_that("MiR_ORA and MiR_Score return enrichment results for the signal pathway", {
  input_data <- get_all_input_data(
    methods = c("MiR_ORA", "MiR_Score"), pathway = toy_pathway,
    mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
    background_MGI = toy_background_MGI
  )

  ora <- MiR_ORA(DEmiR = input_data$data$MiR_ORA, background = unique(input_data$pathway$MiR$miRNA),
                 pathway = input_data$pathway$MiR, pvalueCutoff = 1)
  expect_s3_class(ora, "data.frame")
  expect_true("REPLICATIVE IMMORTALITY" %in% ora$pathway)

  score <- MiR_Score(miRList = input_data$data$MiR_Score, pathway = input_data$pathway$MiR,
                      pvalueType = "std", pvalueCutoff = 1)
  expect_s3_class(score, "data.frame")
  expect_true("REPLICATIVE IMMORTALITY" %in% score$pathway)
})

test_that("MiR_ORA rejects an empty DEmiR vector", {
  expect_error(MiR_ORA(DEmiR = character(0), pathway = toy_pathway), "must be entered")
})
