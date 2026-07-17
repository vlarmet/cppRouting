# aggregate_aux: single legacy vector, single named column out of a
# data.frame, and multi-column data.frame - across get_distance_pair(),
# get_distance_matrix(), normal/CH/CCH graphs and cpp_simplify().

# A strictly linear, tie-free graph: only one possible path between any two
# nodes, so aggregate_aux values can be hand-verified exactly.
line_edges <- function() {
  data.frame(from = c(1, 2, 3), to = c(2, 3, 4), dist = c(4, 6, 2))
}

test_that("legacy single aux vector: aggregate_aux = TRUE returns a plain vector/matrix", {
  e <- line_edges()
  g <- makegraph(e, directed = TRUE, aux = e$dist * 2) # aux = 2x the main weight

  d <- get_distance_pair(g, 1, 4)
  a <- get_distance_pair(g, 1, 4, aggregate_aux = TRUE)
  expect_true(is.numeric(a) && !is.list(a))
  expect_equal(a, 2 * d)

  m <- get_distance_matrix(g, c(1, 2), c(3, 4))
  am <- get_distance_matrix(g, c(1, 2), c(3, 4), aggregate_aux = TRUE)
  expect_true(is.matrix(am))
  expect_equal(unname(am), unname(2 * m))
})

test_that("data.frame aux with a single column resolved returns a length-1 list, not a bare vector/matrix", {
  e <- line_edges()
  g <- makegraph(e, directed = TRUE, aux = data.frame(dist2 = e$dist * 2))

  a <- get_distance_pair(g, 1, 4, aggregate_aux = TRUE)
  expect_s3_class(a, "data.frame")
  expect_named(a, "dist2")
  expect_equal(a$dist2, 2 * get_distance_pair(g, 1, 4))

  am <- get_distance_matrix(g, c(1, 2), c(3, 4), aggregate_aux = TRUE)
  expect_type(am, "list")
  expect_named(am, "dist2")
  expect_true(is.matrix(am$dist2))
})

test_that("multi-column data.frame aux: aggregate_aux = TRUE returns all columns, named", {
  e <- line_edges()
  g <- makegraph(e, directed = TRUE, aux = data.frame(hop = rep(1, nrow(e)), dist2 = e$dist * 2))
  d <- get_distance_pair(g, 1, 4)

  a <- get_distance_pair(g, 1, 4, aggregate_aux = TRUE)
  expect_s3_class(a, "data.frame")
  expect_named(a, c("hop", "dist2"))
  expect_equal(a$hop, 3) # 3 edges on the only path 1->2->3->4
  expect_equal(a$dist2, 2 * d)

  am <- get_distance_matrix(g, c(1, 2), c(3, 4), aggregate_aux = TRUE)
  expect_type(am, "list")
  expect_named(am, c("hop", "dist2"))
})

test_that("character vector aggregate_aux selects a subset of columns", {
  e <- line_edges()
  g <- makegraph(e, directed = TRUE, aux = data.frame(hop = rep(1, nrow(e)), dist2 = e$dist * 2, dist3 = e$dist * 3))

  a <- get_distance_pair(g, 1, 4, aggregate_aux = c("dist2", "dist3"))
  expect_named(a, c("dist2", "dist3"))
  expect_false("hop" %in% names(a))

  a1 <- get_distance_pair(g, 1, 4, aggregate_aux = "hop")
  expect_s3_class(a1, "data.frame")
  expect_named(a1, "hop")
})

test_that("aggregate_aux errors on unknown column names and on a graph without aux", {
  e <- line_edges()
  g_no_aux <- makegraph(e, directed = TRUE)
  expect_error(get_distance_pair(g_no_aux, 1, 4, aggregate_aux = TRUE))

  g <- makegraph(e, directed = TRUE, aux = data.frame(hop = rep(1, nrow(e))))
  expect_error(get_distance_pair(g, 1, 4, aggregate_aux = "not_a_column"))
})

test_that("negative aux values trigger a warning at makegraph() time", {
  e <- line_edges()
  expect_warning(makegraph(e, directed = TRUE, aux = c(-1, 2, 3)))
  expect_warning(makegraph(e, directed = TRUE, aux = data.frame(a = c(-1, 2, 3))))
})

test_that("multi-column aggregate_aux agrees between normal graph and CH (no ties)", {
  g <- grid_graph(6, 6, directed = TRUE, with_aux = TRUE)
  ch <- cpp_contract(g, silent = TRUE)

  set.seed(11)
  nodes <- g$dict$ref
  from <- sample(nodes, 20, replace = TRUE)
  to <- sample(nodes, 20, replace = TRUE)

  a_normal <- get_distance_pair(g, from, to, aggregate_aux = TRUE)
  a_ch <- get_distance_pair(ch, from, to, aggregate_aux = TRUE)
  expect_equal(a_normal$hop, a_ch$hop, tolerance = 1e-8)
  expect_equal(a_normal$dist2, a_ch$dist2, tolerance = 1e-6)

  m_normal <- get_distance_matrix(g, from[1:6], to[1:6], aggregate_aux = TRUE)
  m_ch <- get_distance_matrix(ch, from[1:6], to[1:6], aggregate_aux = TRUE)
  expect_equal(unname(m_normal$hop), unname(m_ch$hop), tolerance = 1e-8)
  expect_equal(unname(m_normal$dist2), unname(m_ch$dist2), tolerance = 1e-6)
})

test_that("aggregate_aux on a tied-cost graph: the aux column proportional to the main weight always matches exactly, even when paths differ", {
  # toy_graph has genuine cost ties (e.g. 1->4 via 1-2-4 or 1-3-4, both cost
  # 5) - documented, expected behavior (see get_distance_pair()'s @details).
  # dist2 = 2 * dist is linear in the optimized weight, so its sum along any
  # optimal path always equals 2 * distance, regardless of which tied path a
  # given algorithm/graph representation happens to pick.
  g <- toy_graph(aux = toy_aux())
  ch <- cpp_contract(g, silent = TRUE)

  nodes <- g$dict$ref
  from <- rep(nodes, times = length(nodes))
  to <- rep(nodes, each = length(nodes))

  d_normal <- get_distance_pair(g, from, to)
  a_normal <- get_distance_pair(g, from, to, aggregate_aux = TRUE)
  a_ch <- get_distance_pair(ch, from, to, aggregate_aux = TRUE)

  ok <- !is.na(d_normal)
  expect_equal(a_normal$dist2[ok], 2 * d_normal[ok], tolerance = 1e-8)
  expect_equal(a_ch$dist2[ok], 2 * d_normal[ok], tolerance = 1e-8)
})

test_that("cpp_simplify sums multi-column aux correctly across merged edge chains", {
  # a-b-c-d is a pure chain (b, c are non-intersection nodes and get
  # contracted away by cpp_simplify); the merged edge's aux must equal the
  # sum over the 3 original edges.
  e <- data.frame(from = c("a", "b", "c"), to = c("b", "c", "d"), dist = c(1, 2, 3))
  aux <- data.frame(hop = c(1, 1, 1), dist2 = e$dist * 2)
  g <- makegraph(e, directed = TRUE, aux = aux)

  simplified <- cpp_simplify(g, iterate = TRUE)
  expect_equal(nrow(simplified$data), 1) # a-b-c-d collapsed to a single a->d edge

  a <- get_distance_pair(simplified, "a", "d", aggregate_aux = TRUE)
  expect_equal(a$hop, 3)
  expect_equal(a$dist2, sum(aux$dist2))
})
