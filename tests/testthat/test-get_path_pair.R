test_that("keep filters the returned path down to the requested nodes (endpoints included only if listed)", {
  g <- grid_graph(4, 4, directed = TRUE)
  full <- get_path_pair(g, "1", "16")[[1]]
  expect_true("9" %in% full) # sanity: node 9 is really on the unrestricted path

  kept <- get_path_pair(g, "1", "16", keep = "9")[[1]]
  expect_equal(kept, "9") # endpoints are not force-included if absent from `keep`

  kept_endpoints <- get_path_pair(g, "1", "16", keep = c("1", "16"))[[1]]
  expect_equal(kept_endpoints, c("1", "16"))
})

test_that("long = TRUE returns a long data.frame with one row per path node, and drops unreachable pairs entirely", {
  e <- data.frame(from = c(1, 2, 4), to = c(2, 3, 5), dist = c(1, 1, 1))
  g <- makegraph(e, directed = TRUE)

  res <- get_path_pair(g, c("1", "1"), c("3", "5"), long = TRUE)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("from", "to", "node"))
  expect_true(all(res$from == "1" & res$to == "3")) # the 1->5 (unreachable) pair contributes 0 rows
  expect_equal(res$node, c("1", "2", "3"))
})

test_that("unreachable pairs return an empty path (character(0)), not an error", {
  e <- data.frame(from = c(1, 2, 4), to = c(2, 3, 5), dist = c(1, 1, 1))
  g <- makegraph(e, directed = TRUE)
  p <- get_path_pair(g, "1", "5")[[1]]
  expect_length(p, 0)
})

test_that("results are named 'from_to' in list form", {
  g <- toy_graph()
  res <- get_path_pair(g, c("1", "2"), c("6", "5"))
  expect_equal(names(res), c("1_6", "2_5"))
})
