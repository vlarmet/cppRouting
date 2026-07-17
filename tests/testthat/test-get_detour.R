# get_detour()'s documented condition: a returned node n satisfies
# SP(o,n) + SP(n,d) < SP(o,d) + extra, with o = from, d = to.

test_that("every returned detour node satisfies the documented distance condition", {
  g <- grid_graph(5, 5, directed = TRUE)
  from <- "1"; to <- "25"; extra <- 5

  base_dist <- get_distance_pair(g, from, to)
  detour_nodes <- get_detour(g, from = from, to = to, extra = extra)[[paste0(from, "_", to)]]
  expect_gt(length(detour_nodes), 0)

  d_from <- get_distance_pair(g, rep(from, length(detour_nodes)), detour_nodes)
  d_to <- get_distance_pair(g, detour_nodes, rep(to, length(detour_nodes)))
  expect_true(all(d_from + d_to < base_dist + extra + 1e-8))
})

test_that("a larger extra budget returns at least as many nodes", {
  g <- grid_graph(5, 5, directed = TRUE)
  small <- get_detour(g, from = "1", to = "25", extra = 2)[["1_25"]]
  large <- get_detour(g, from = "1", to = "25", extra = 10)[["1_25"]]
  expect_true(all(small %in% large))
})

test_that("nodes on the shortest path itself are always included (extra detour cost 0)", {
  g <- grid_graph(5, 5, directed = TRUE)
  path <- get_path_pair(g, "1", "25")[["1_25"]]
  detour_nodes <- get_detour(g, from = "1", to = "25", extra = 1)[["1_25"]]
  expect_true(all(path %in% detour_nodes))
})
