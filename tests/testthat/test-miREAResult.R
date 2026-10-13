test_that("miREAResult constructor and accessors work", {
  res <- miREAResult(
    result = list(TG_ORA = data.frame(pathway = "PW1", padj = 0.01)),
    time = data.frame(method = "TG_ORA", time = 0.1),
    valid_methods = "TG_ORA",
    parameter = list(n_pathway = 1)
  )

  expect_s4_class(res, "miREAResult")
  expect_equal(res$valid_methods, "TG_ORA")
  expect_equal(res[["valid_methods"]], "TG_ORA")
  expect_equal(res$time$time, 0.1)
  expect_equal(res[["parameter"]]$n_pathway, 1)
  expect_output(print(res), "miREAResult object")

  df <- as.data.frame(res)
  expect_equal(df$method, "TG_ORA")
  expect_equal(nrow(df), 1)
})
