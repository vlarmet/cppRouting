test_that("default weights match the original graph cost", {
  g <- toy_graph()
  cch <- cpp_contract(g, customizable = TRUE, silent = TRUE)
  metric <- cpp_customize(cch)
  expect_s3_class(metric, "cppRouting_cch_metric")

  d_metric <- get_distance_pair(metric, "1", "6")
  d_normal <- get_distance_pair(g, "1", "6")
  expect_equal(d_metric, d_normal, tolerance = 1e-8)
})

test_that("custom weights change the result as expected", {
  g <- toy_graph()
  cch <- cpp_contract(g, customizable = TRUE, silent = TRUE)

  original_weights <- g$data$dist
  inflated_weights <- original_weights * 10 # scaling every edge inflates every path by the same factor
  metric_default <- cpp_customize(cch)
  metric_inflated <- cpp_customize(cch, weights = inflated_weights)

  d_default <- get_distance_pair(metric_default, "1", "6")
  d_inflated <- get_distance_pair(metric_inflated, "1", "6")
  expect_equal(d_inflated, d_default * 10, tolerance = 1e-6)
})

test_that("weights are validated: wrong length, NA and negative are rejected", {
  g <- toy_graph()
  cch <- cpp_contract(g, customizable = TRUE, silent = TRUE)
  n <- nrow(g$data)

  expect_error(cpp_customize(cch, weights = rep(1, n - 1)))
  expect_error(cpp_customize(cch, weights = c(NA, rep(1, n - 1))))
  expect_error(cpp_customize(cch, weights = c(-1, rep(1, n - 1))))
})

test_that("only a cpp_contract(customizable = TRUE) object is accepted", {
  g <- toy_graph()
  ch <- cpp_contract(g, silent = TRUE) # classical, not customizable
  expect_error(cpp_customize(ch))
  expect_error(cpp_customize(g))
})

test_that("querying an un-customized CCH topology errors with a clear message", {
  g <- toy_graph()
  cch <- cpp_contract(g, customizable = TRUE, silent = TRUE)
  expect_error(get_distance_pair(cch, "1", "6"), "cpp_customize")
})
