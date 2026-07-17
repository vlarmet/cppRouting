test_that("get_aon returns a data.frame with from/to/cost/flow, non-negative flow", {
  g <- toy_graph()
  res <- get_aon(g, from = "1", to = "6", demand = 10)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("from", "to", "cost", "flow"))
  expect_true(all(res$flow >= 0))
})

test_that("flow is conserved: all demand leaves the origin and arrives at the destination", {
  # robust to tie-breaking on the shortest path itself (see aggregate_aux
  # tests for why toy_graph has genuine ties): whichever path is used, the
  # total flow crossing out of the origin / into the destination must equal
  # the demand.
  g <- toy_graph()
  res <- get_aon(g, from = "1", to = "6", demand = 10)
  expect_equal(sum(res$flow[res$from == "1"]), 10)
  expect_equal(sum(res$flow[res$to == "6"]), 10)
})

test_that("get_aon works on normal, CH and CCH graphs and gives consistent total flow", {
  g <- grid_graph(5, 5, directed = TRUE)
  ch <- cpp_contract(g, silent = TRUE)
  cch <- cpp_customize(cpp_contract(g, customizable = TRUE, silent = TRUE))

  from <- c("1", "2"); to <- c("25", "24"); demand <- c(5, 7)

  res_normal <- get_aon(g, from, to, demand)
  res_ch <- get_aon(ch, from, to, demand)
  res_cch <- get_aon(cch, from, to, demand)

  expect_equal(sum(res_normal$flow), sum(res_ch$flow))
  expect_equal(sum(res_normal$flow), sum(res_cch$flow))
})

test_that("invalid algorithm values are rejected", {
  g <- toy_graph()
  expect_error(get_aon(g, "1", "6", 10, algorithm = "not_an_algorithm"))
})
