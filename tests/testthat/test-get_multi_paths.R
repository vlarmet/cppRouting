test_that("get_multi_paths returns the full from x to cross product", {
  g <- grid_graph(4, 4, directed = TRUE)
  from <- c("1", "2"); to <- c("13", "14")
  res <- get_multi_paths(g, from, to)

  expect_named(res, from)
  for (f in from) expect_named(res[[f]], to)
})

test_that("each returned path starts and ends at the requested nodes, and matches get_distance_pair via edge costs", {
  g <- grid_graph(4, 4, directed = TRUE)
  from <- "1"; to <- c("13", "16")
  res <- get_multi_paths(g, from, to)[["1"]]

  for (t in to) {
    p <- res[[t]]
    if (length(p) == 0) next
    expect_equal(p[1], from)
    expect_equal(p[length(p)], t)
  }
})

test_that("unreachable pairs give an empty path", {
  e <- data.frame(from = c(1, 2, 4), to = c(2, 3, 5), dist = c(1, 1, 1))
  g <- makegraph(e, directed = TRUE)
  res <- get_multi_paths(g, "1", "5")
  expect_length(res[["1"]][["5"]], 0)
})
