# Cross-graph-type agreement: normal graph, CH (classical contraction), CCH
# (customizable contraction) and a simplified graph must all agree on
# shortest-path distances for the same origin/destination pairs. This is the
# regression net that would have caught most of the correctness bugs fixed
# during CCH integration (get_path_pair on CCH, min_degree_order, the
# distance_matrix_cch bucket rewrite, etc.).
#
# Continuous random edge weights (see helper-graphs.R) make exact ties
# between distinct shortest paths astronomically unlikely, so distances are
# expected to match exactly (up to floating point noise) - not just "close".

# Cheapest edge per (from, to), matching what makegraph()/the routing
# engines do when parallel edges exist (used to verify a returned path's own
# cost against the distance reported for the same pair).
edge_lookup <- function(edges) {
  key <- paste(edges$from, edges$to)
  ord <- order(edges$dist)
  key_sorted <- key[ord]
  keep <- !duplicated(key_sorted)
  stats::setNames(edges$dist[ord][keep], key_sorted[keep])
}

path_cost <- function(path, lookup) {
  if (length(path) < 2 || any(is.na(path))) return(NA_real_)
  keys <- paste(utils::head(path, -1), utils::tail(path, -1))
  sum(unname(lookup[keys]))
}

test_that("normal, CH and CCH graphs agree on distances", {
  g <- grid_graph(6, 6, directed = TRUE)
  ch <- cpp_contract(g, silent = TRUE)
  cch <- cpp_customize(cpp_contract(g, customizable = TRUE, silent = TRUE))

  set.seed(2)
  nodes <- g$dict$ref
  from <- sample(nodes, 30, replace = TRUE)
  to <- sample(nodes, 30, replace = TRUE)

  d_normal <- get_distance_pair(g, from, to)
  d_ch <- get_distance_pair(ch, from, to)
  d_cch <- get_distance_pair(cch, from, to)

  expect_equal(d_normal, d_ch, tolerance = 1e-8)
  expect_equal(d_normal, d_cch, tolerance = 1e-8)
})

test_that("the 4 unconstracted-graph algorithms agree on distances", {
  g <- grid_graph(6, 6, directed = TRUE, with_coords = TRUE)
  set.seed(3)
  nodes <- g$dict$ref
  from <- sample(nodes, 25, replace = TRUE)
  to <- sample(nodes, 25, replace = TRUE)

  d_dijkstra <- get_distance_pair(g, from, to, algorithm = "Dijkstra")
  d_bi <- get_distance_pair(g, from, to, algorithm = "bi")
  # constant = 1 is intentional here: edge weights are >= 1 (see grid_edges())
  # and grid steps are unit distance, so the euclidean heuristic is
  # admissible; suppress the "are you sure?" nudge get_distance_pair() gives
  # by default for constant == 1.
  d_astar <- suppressWarnings(get_distance_pair(g, from, to, algorithm = "A*"))
  d_nba <- suppressWarnings(get_distance_pair(g, from, to, algorithm = "NBA"))

  expect_equal(d_dijkstra, d_bi, tolerance = 1e-8)
  expect_equal(d_dijkstra, d_astar, tolerance = 1e-8)
  expect_equal(d_dijkstra, d_nba, tolerance = 1e-8)
})

test_that("get_distance_matrix agrees across normal, mch and phast", {
  g <- grid_graph(6, 6, directed = TRUE)
  ch <- cpp_contract(g, silent = TRUE)

  set.seed(4)
  nodes <- g$dict$ref
  from <- sample(nodes, 8)
  to <- sample(nodes, 10)

  m_normal <- get_distance_matrix(g, from, to)
  m_mch <- get_distance_matrix(ch, from, to, algorithm = "mch")
  m_phast <- get_distance_matrix(ch, from, to, algorithm = "phast")

  expect_equal(unname(m_normal), unname(m_mch), tolerance = 1e-8)
  expect_equal(unname(m_normal), unname(m_phast), tolerance = 1e-8)
})

test_that("simplified graph preserves distances between kept nodes", {
  g <- grid_graph(6, 6, directed = TRUE)
  set.seed(5)
  nodes <- g$dict$ref
  from <- sample(nodes, 15, replace = TRUE)
  to <- sample(nodes, 15, replace = TRUE)
  keep_nodes <- unique(c(from, to))

  simplified <- cpp_simplify(g, iterate = TRUE, keep = keep_nodes)
  expect_s3_class(simplified, "cppRouting_simplified")
  expect_s3_class(simplified, "cppRouting_graph")

  d_normal <- get_distance_pair(g, from, to)
  d_simplified <- get_distance_pair(simplified, from, to)
  expect_equal(d_normal, d_simplified, tolerance = 1e-8)
})

test_that("get_path_pair returns paths whose own cost matches get_distance_pair, on all graph types", {
  g <- grid_graph(6, 6, directed = TRUE)
  edges <- to_df(g)
  colnames(edges) <- c("from", "to", "dist")
  lookup <- edge_lookup(edges)

  ch <- cpp_contract(g, silent = TRUE)
  cch <- cpp_customize(cpp_contract(g, customizable = TRUE, silent = TRUE))

  set.seed(6)
  nodes <- g$dict$ref
  from <- sample(nodes, 20, replace = TRUE)
  to <- sample(nodes, 20, replace = TRUE)

  for (graph in list(normal = g, ch = ch, cch = cch)) {
    d <- get_distance_pair(graph, from, to)
    paths <- get_path_pair(graph, from, to)
    costs <- vapply(unname(paths), path_cost, numeric(1), lookup = lookup)
    expect_equal(unname(d), costs, tolerance = 1e-6)
  }
})

test_that("get_path_pair paths start and end at the requested nodes", {
  g <- grid_graph(6, 6, directed = TRUE)
  set.seed(7)
  nodes <- g$dict$ref
  from <- sample(nodes, 10, replace = TRUE)
  to <- sample(nodes, 10, replace = TRUE)

  paths <- get_path_pair(g, from, to)
  for (i in seq_along(paths)) {
    p <- paths[[i]]
    if (length(p) >= 1 && !anyNA(p)) {
      expect_equal(as.character(p[1]), as.character(from[i]))
      expect_equal(as.character(p[length(p)]), as.character(to[i]))
    }
  }
})

test_that("self-pairs (from == to) have distance 0 and a single-node path", {
  g <- grid_graph(6, 6, directed = TRUE)
  ch <- cpp_contract(g, silent = TRUE)
  cch <- cpp_customize(cpp_contract(g, customizable = TRUE, silent = TRUE))
  node <- g$dict$ref[1]

  expect_equal(get_distance_pair(g, node, node), 0)
  expect_equal(get_distance_pair(ch, node, node), 0)
  expect_equal(get_distance_pair(cch, node, node), 0)

  p <- get_path_pair(g, node, node)[[1]]
  expect_equal(as.character(p), as.character(node))
})

test_that("parallel edges: the cheapest one is used, consistently across graph types", {
  edges <- with_parallel_edge(grid_edges(6, 6), from = 1, to = 2, extra_cost = 50)
  g <- makegraph(edges, directed = TRUE)
  ch <- cpp_contract(g, silent = TRUE)

  d_normal <- get_distance_pair(g, "1", "2")
  d_ch <- get_distance_pair(ch, "1", "2")
  cheapest_direct <- min(edges$weight[edges$from == 1 & edges$to == 2])

  expect_lte(d_normal, cheapest_direct + 1e-8)
  expect_equal(d_normal, d_ch, tolerance = 1e-8)
})
