test_that("matrix has the right dimensions and dimnames", {
  g <- grid_graph(5, 5, directed = TRUE)
  from <- as.character(c(1, 2, 3))
  to <- as.character(c(10, 11, 12, 13))

  m <- get_distance_matrix(g, from, to)
  expect_equal(dim(m), c(length(from), length(to)))
  expect_equal(rownames(m), from)
  expect_equal(colnames(m), to)
})

test_that("matrix is consistent regardless of which side (from/to) is larger", {
  # get_distance_matrix internally roots the search at whichever of from/to
  # is smaller and transposes if needed - both orientations must agree.
  # This requires an undirected (symmetric) graph: on a directed one,
  # distance(A->B) and distance(B->A) are generally different values, not
  # transposes of each other, so grid_graph()'s directed-only-right/down
  # topology would not be a valid check here.
  g <- grid_graph(6, 6, directed = FALSE)
  small <- as.character(1:3)
  large <- as.character(10:20)

  m1 <- get_distance_matrix(g, small, large) # from smaller than to
  m2 <- get_distance_matrix(g, large, small) # from larger than to
  expect_equal(unname(m1), unname(t(m2)), tolerance = 1e-8)
})

test_that("CH matrix orientation is also consistent (mch and phast)", {
  g <- grid_graph(6, 6, directed = FALSE)
  ch <- cpp_contract(g, silent = TRUE)
  small <- as.character(1:3)
  large <- as.character(10:20)

  for (algo in c("mch", "phast")) {
    m1 <- get_distance_matrix(ch, small, large, algorithm = algo)
    m2 <- get_distance_matrix(ch, large, small, algorithm = algo)
    expect_equal(unname(m1), unname(t(m2)), tolerance = 1e-8)
  }
})

test_that("unreachable pairs are NA in the matrix", {
  e <- data.frame(from = c(1, 2, 4), to = c(2, 3, 5), dist = c(1, 1, 1))
  g <- makegraph(e, directed = TRUE)
  m <- get_distance_matrix(g, from = c("1"), to = c("3", "5"))
  expect_equal(unname(m[1, 1]), 2)
  expect_true(is.na(m[1, 2]))
})

test_that("legacy plain-list graphs (no class) are still accepted", {
  # cppRouting < 3.0 and some downstream packages (e.g. rTARDIS) pass plain
  # lists without the cppRouting_graph/cppRouting_contracted class.
  g <- grid_graph(5, 5, directed = TRUE)
  legacy <- as.list(g)
  class(legacy) <- "list"
  m <- get_distance_matrix(legacy, from = c("1"), to = c("10", "11"))
  expect_equal(unname(m), unname(get_distance_matrix(g, c("1"), c("10", "11"))))

  ch <- cpp_contract(g, silent = TRUE)
  legacy_ch <- as.list(ch)
  class(legacy_ch) <- "list"
  m2 <- get_distance_matrix(legacy_ch, from = c("1"), to = c("10", "11"))
  expect_equal(unname(m2), unname(get_distance_matrix(ch, c("1"), c("10", "11"))))
})
