data(toy_pathway, toy_gene_DEdata, toy_mir_DEdata, toy_gene_mat, toy_mir_mat,
     toy_background_MGI, toy_background_GGI, package = "miREA")

test_that("miREA() runs all nine methods and returns a result per method", {
  input_data <- suppressWarnings(get_all_input_data(
    methods = default.method_list, pathway = toy_pathway,
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
    methods = default.method_list, input_data = input_data, background = background,
    minSize = 2, maxSize = 5000,
    pvalueType = list(TG = "std", MiR = "std", Edge = "neg"),
    pAdjMethod = "BH", pvalueCutoff = 1, iter = 100, ncores = 1
  ))

  expect_setequal(result$valid_methods, default.method_list)
  for (m in default.method_list) {
    expect_s3_class(result$result[[m]], "data.frame")
    expect_equal(nrow(result$result[[m]]), 3)
  }
})

test_that("miREA() requires a background when none is supplied", {
  input_data <- suppressWarnings(get_all_input_data(
    methods = "TG_ORA", pathway = toy_pathway,
    mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
    background_MGI = toy_background_MGI
  ))
  expect_error(
    miREA(methods = "TG_ORA", input_data = input_data, background = NULL),
    "background must be supplied"
  )
})
