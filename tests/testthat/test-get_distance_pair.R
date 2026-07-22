test_that("from and to must have the same length", {
  g <- toy_graph()
  expect_error(get_distance_pair(g, from = c(1, 2), to = 1))
})

test_that("NAs in from/to are rejected, unknown nodes are rejected", {
  g <- toy_graph()
  expect_error(get_distance_pair(g, from = NA, to = 1))
  expect_error(get_distance_pair(g, from = "not_a_node", to = 1))
})

test_that("unreachable pairs return NA, not an error", {
  # two disconnected components: 1-2-3 and 4-5, no bridge between them
  e <- data.frame(from = c(1, 2, 4), to = c(2, 3, 5), dist = c(1, 1, 1))
  g <- makegraph(e, directed = TRUE)
  d <- get_distance_pair(g, from = c(1, 1), to = c(3, 5))
  expect_equal(d[1], 1 + 1)
  expect_true(is.na(d[2]))
})

test_that("vectorized calls match calling one pair at a time", {
  g <- toy_graph()
  nodes <- g$dict$ref
  set.seed(21)
  from <- sample(nodes, 8, replace = TRUE)
  to <- sample(nodes, 8, replace = TRUE)

  batched <- get_distance_pair(g, from, to)
  single <- mapply(function(f, t) get_distance_pair(g, f, t), from, to)
  expect_equal(unname(batched), unname(single))
})

test_that("A*/NBA warn when constant == 1 with coordinates but still return a value", {
  g <- grid_graph(4, 4, directed = TRUE, with_coords = TRUE)
  expect_warning(d <- get_distance_pair(g, "1", "16", algorithm = "A*"))
  expect_true(is.numeric(d) && !is.na(d))
})
