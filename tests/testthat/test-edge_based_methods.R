data(toy_pathway, toy_gene_DEdata, toy_mir_DEdata, toy_gene_mat, toy_mir_mat,
     toy_background_MGI, toy_background_GGI, package = "miREA")

test_that("Edge_ORA, Edge_Score1D, Edge_Score2D return enrichment results for the signal pathway", {
  input_data <- get_all_input_data(
    methods = c("Edge_ORA", "Edge_Score1D", "Edge_Score2D"), pathway = toy_pathway,
    mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
    background_MGI = toy_background_MGI, gene_mat = toy_gene_mat, mir_mat = toy_mir_mat,
    scoreFun = "rank"
  )
  background <- unique(paste0(toy_background_MGI$miRNA, ":", toy_background_MGI$gene))

  ora <- Edge_ORA(DEMGI = input_data$data$Edge_ORA, background = background,
                   pathway = input_data$pathway$Edge, pvalueCutoff = 1)
  expect_s3_class(ora, "data.frame")
  expect_true("REPLICATIVE IMMORTALITY" %in% ora$pathway)

  score1d <- Edge_Score1D(MGIList = input_data$data$Edge_Score1D, pathway = input_data$pathway$Edge,
                           pvalueType = "neg", pvalueCutoff = 1)
  expect_s3_class(score1d, "data.frame")
  expect_true("REPLICATIVE IMMORTALITY" %in% score1d$pathway)

  score2d <- Edge_Score2D(MGI2D = input_data$data$Edge_Score2D, pathway = input_data$pathway$Edge,
                           background = background, iter = 100, ncores = 1, pvalueCutoff = 1)
  expect_s3_class(score2d, "data.frame")
  expect_true("REPLICATIVE IMMORTALITY" %in% score2d$pathway)
})

test_that("Edge_ORA rejects an empty DEMGI vector", {
  expect_error(Edge_ORA(DEMGI = character(0), pathway = toy_pathway), "DEMGI is empty")
})

test_that("Edge_Topology returns enrichment results for the signal pathway", {
  input_data <- get_all_input_data(
    methods = "Edge_Topology", pathway = toy_pathway,
    mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
    background_MGI = toy_background_MGI, background_GGI = toy_background_GGI,
    GGI_source = "Omnipath", gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
  )
  background <- unique(paste0(toy_background_MGI$miRNA, ":", toy_background_MGI$gene))

  res <- Edge_Topology(
    MGI = input_data$data$Edge_Topology, background = background,
    pathway = input_data$pathway$Edge, pathway_GGI = input_data$pathway_GGI,
    dict_MGI = input_data$dict_MGI, pvalueType = "neg",
    minSize = 2, maxSize = 5000, iter = 100, ncores = 1, pvalueCutoff = 1
  )
  expect_s3_class(res, "data.frame")
  expect_true("REPLICATIVE IMMORTALITY" %in% res$pathway)
})

test_that("Edge_Network returns enrichment results for the signal pathway", {
  input_data <- get_all_input_data(
    methods = "Edge_Network", pathway = toy_pathway,
    mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
    background_MGI = toy_background_MGI, background_GGI = toy_background_GGI,
    GGI_source = "Omnipath", gene_mat = toy_gene_mat, mir_mat = toy_mir_mat
  )

  res <- Edge_Network(
    DEMGI = input_data$data$Edge_Network, pathway_gene = input_data$pathway$TG,
    pathway_MGI = input_data$pathway$Edge, pathway_GGI = input_data$pathway_GGI,
    MGI = input_data$MGI, GGI = input_data$GGI,
    minSize = 2, maxSize = 5000, iter = 100, ncores = 1, pvalueCutoff = 1
  )
  expect_s3_class(res, "data.frame")
  expect_true("REPLICATIVE IMMORTALITY" %in% res$pathway)
})
