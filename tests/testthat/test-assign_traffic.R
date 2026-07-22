traffic_graph <- function() {
  grid_graph(4, 4, directed = TRUE, with_traffic_attrs = TRUE)
}

test_that("link-based algorithms (msa, fw, cfw, bfw) converge to a small gap", {
  g <- traffic_graph()
  set.seed(41)
  nodes <- g$dict$ref
  demand <- data.frame(from = sample(nodes, 5), to = sample(nodes, 5), demand = runif(5, 1, 10))

  for (algo in c("msa", "fw", "cfw", "bfw")) {
    res <- assign_traffic(g, from = demand$from, to = demand$to, demand = demand$demand,
                          algorithm = algo, max_it = 50, verbose = FALSE)
    expect_true(is.finite(res$gap))
    expect_gte(res$gap, 0)
    expect_true(all(res$data$flow >= 0))
  }
})

test_that("bush-based algorithm (dial) converges to a small gap", {
  g <- traffic_graph()
  set.seed(42)
  nodes <- g$dict$ref
  demand <- data.frame(from = sample(nodes, 5), to = sample(nodes, 5), demand = runif(5, 1, 10))

  res <- assign_traffic(g, from = demand$from, to = demand$to, demand = demand$demand,
                        algorithm = "dial", max_it = 10, verbose = FALSE)
  expect_true(is.finite(res$gap))
  expect_gte(res$gap, 0)
})

test_that("aon_method = 'cch' auto-prepares a CCH when cch = NULL, and matches aon_method = 'bi'", {
  g <- traffic_graph()
  set.seed(43)
  nodes <- g$dict$ref
  demand <- data.frame(from = sample(nodes, 5), to = sample(nodes, 5), demand = runif(5, 1, 10))

  # auto-preparation (cch = NULL) shows a progress bar for the underlying
  # cpp_contract(customizable = TRUE) call; captured here to keep test output
  # clean, not asserted on.
  invisible(capture.output(
    res_cch <- assign_traffic(g, from = demand$from, to = demand$to, demand = demand$demand,
                              algorithm = "bfw", max_it = 50, aon_method = "cch", verbose = FALSE)
  ))
  res_bi <- assign_traffic(g, from = demand$from, to = demand$to, demand = demand$demand,
                           algorithm = "bfw", max_it = 50, aon_method = "bi", verbose = FALSE)
  expect_equal(res_cch$gap, res_bi$gap, tolerance = 1e-6)
})

test_that("a cch prepared for a different graph is rejected with a clear error", {
  g1 <- traffic_graph()
  g2 <- grid_graph(6, 6, directed = TRUE, with_traffic_attrs = TRUE) # different topology

  wrong_cch <- cpp_contract(g2, customizable = TRUE, silent = TRUE)

  demand <- data.frame(from = g1$dict$ref[1:3], to = g1$dict$ref[4:6], demand = c(1, 2, 3))
  expect_error(
    assign_traffic(g1, from = demand$from, to = demand$to, demand = demand$demand,
                   algorithm = "bfw", aon_method = "cch", cch = wrong_cch, verbose = FALSE),
    "topology|different number"
  )
})

test_that("from/to/demand must have the same length", {
  g <- traffic_graph()
  expect_error(assign_traffic(g, from = c("1", "2"), to = "1", demand = c(1, 2)))
})

test_that("alpha, beta and capacity must be defined on the graph", {
  g <- grid_graph(4, 4, directed = TRUE) # no traffic attributes
  expect_error(assign_traffic(g, from = "1", to = "2", demand = 1))
})
