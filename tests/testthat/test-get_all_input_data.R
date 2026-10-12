data(toy_pathway, toy_gene_DEdata, toy_mir_DEdata, toy_gene_mat, toy_mir_mat,
     toy_background_MGI, toy_background_GGI, package = "miREA")

test_that("get_all_input_data() builds all expected components for the full method list", {
  input_data <- get_all_input_data(
    methods = default.method_list,
    pathway = toy_pathway,
    mir_DEdata = toy_mir_DEdata,
    gene_DEdata = toy_gene_DEdata,
    background_MGI = toy_background_MGI,
    background_GGI = toy_background_GGI,
    GGI_source = "Omnipath",
    gene_mat = toy_gene_mat,
    mir_mat = toy_mir_mat,
    scoreFun = "rank"
  )

  expect_type(input_data, "list")
  expect_true(all(c("data", "pathway", "pathway_GGI", "GGI", "MGI", "dict_MGI") %in% names(input_data)))
  expect_true(all(default.method_list %in% names(input_data$data)))
  expect_s3_class(input_data$dict_MGI, "data.frame")
  expect_true(nrow(input_data$dict_MGI) > 0)
})

test_that("get_all_input_data() requires background_MGI", {
  expect_error(
    get_all_input_data(
      methods = "TG_ORA", pathway = toy_pathway,
      mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata
    ),
    "background_MGI must be supplied"
  )
})

test_that("get_all_input_data() works for a single method without the full pipeline", {
  input_data <- get_all_input_data(
    methods = "TG_ORA", pathway = toy_pathway,
    mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
    background_MGI = toy_background_MGI
  )
  expect_true(is.character(input_data$data$TG_ORA))
  expect_true(length(input_data$data$TG_ORA) > 0)
})
