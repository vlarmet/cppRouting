test_that("to_igraph() builds a directed igraph with the right vertex/edge counts and names", {
  skip_if_not_installed("igraph")
  g <- toy_graph()
  ig <- to_igraph(g)
  expect_s3_class(ig, "igraph")
  expect_true(igraph::is_directed(ig))
  expect_equal(igraph::vcount(ig), g$nbnode)
  expect_equal(igraph::ecount(ig), nrow(g$data))
  expect_setequal(igraph::V(ig)$name, g$dict$ref)
})

test_that("to_igraph()'s edge weight matches the graph's dist, usable by igraph's own shortest-path functions", {
  skip_if_not_installed("igraph")
  e <- toy_edges()
  g <- makegraph(e, directed = TRUE)
  ig <- to_igraph(g)
  expect_true("weight" %in% igraph::edge_attr_names(ig))
  expect_false("dist" %in% igraph::edge_attr_names(ig)) # renamed, not duplicated

  d_cpp <- get_distance_pair(g, from = "1", to = "6")
  d_ig <- igraph::distances(ig, v = "1", to = "6", mode = "out")[1, 1]
  expect_equal(d_cpp, d_ig, tolerance = 1e-8)
})

test_that("to_igraph() attaches x/y vertex attributes when the graph has coordinates, none otherwise", {
  skip_if_not_installed("igraph")
  e <- toy_edges()
  coords <- data.frame(node = unique(c(e$from, e$to)))
  coords$X <- seq_len(nrow(coords))
  coords$Y <- rev(seq_len(nrow(coords)))
  g <- makegraph(e, directed = TRUE, coords = coords)

  ig <- to_igraph(g)
  expect_true(all(c("x", "y") %in% igraph::vertex_attr_names(ig)))
  expect_equal(igraph::V(ig)$x[match(g$dict$ref, igraph::V(ig)$name)], g$coords[[2]])

  ig_nocoords <- to_igraph(toy_graph()) # toy_graph() has no coords
  expect_false("x" %in% igraph::vertex_attr_names(ig_nocoords))
})

test_that("to_igraph(aux/capacity/alpha/beta = TRUE) carries them over as edge attributes", {
  skip_if_not_installed("igraph")
  e <- toy_edges()
  g <- makegraph(e, directed = TRUE, aux = toy_aux(), capacity = 10, alpha = 0.15, beta = 4)

  ig <- to_igraph(g, aux = TRUE, capacity = TRUE, alpha = TRUE, beta = TRUE)
  expect_true(all(c("hop", "dist2", "capacity", "alpha", "beta") %in% igraph::edge_attr_names(ig)))

  ig_default <- to_igraph(g) # none requested: none of these attributes present
  expect_false(any(c("hop", "dist2", "capacity", "alpha", "beta") %in% igraph::edge_attr_names(ig_default)))
})

test_that("to_igraph() rejects a non-cppRouting_graph input and requires the igraph package", {
  expect_error(to_igraph(list(a = 1)), "must be an object generated")
})
