data(toy_pathway, toy_gene_DEdata, toy_mir_DEdata, toy_gene_mat, toy_mir_mat,
     toy_background_MGI, toy_background_GGI, package = "miREA")

test_that("plotSummary(), plotHeatmap(), and plotHeatmapSankey() run without error", {
  input_data <- suppressWarnings(getAllInputData(
    methods = default_method_list, pathway = toy_pathway,
    mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
    background_MGI = toy_background_MGI, background_GGI = toy_background_GGI,
    GGI_source = "Omnipath", gene_mat = toy_gene_mat, mir_mat = toy_mir_mat,
    scoreFun = "rank"
  ))
  background <- list(
    TG = unique(toy_pathway$gene),
    MiR = unique(input_data$pathway$MiR$miRNA),
    Edge = unique(paste0(toy_background_MGI$miRNA, ":", toy_background_MGI$gene))
  )
  result <- suppressWarnings(miREA(
    methods = default_method_list, input_data = input_data, background = background,
    minSize = 2, maxSize = 5000,
    pvalueType = list(TG = "std", MiR = "std", Edge = "neg"),
    pAdjMethod = "BH", pvalueCutoff = 1, iter = 20, ncores = 1
  ))

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  expect_no_error(plotSummary(result, penrichCutoff = 1))
  expect_no_error(plotHeatmap(method = "Edge_ORA", result = result, input_data = input_data,
                                penrichCutoff = 1))
  expect_no_error(plotHeatmapSankey(method = "Edge_ORA", result = result, input_data = input_data,
                                       penrichCutoff = 1))
})
