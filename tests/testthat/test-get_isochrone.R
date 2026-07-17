test_that("every returned node is actually within the requested time limit", {
  g <- grid_graph(5, 5, directed = TRUE)
  iso <- get_isochrone(g, from = "1", lim = c(3, 6, 9))[["1"]]

  for (lim_name in names(iso)) {
    nodes <- iso[[lim_name]]
    if (length(nodes) == 0) next
    d <- get_distance_pair(g, rep("1", length(nodes)), nodes)
    expect_true(all(d <= as.numeric(lim_name) + 1e-8))
  }
})

test_that("buckets are cumulative by default (nested) and grow monotonically", {
  g <- grid_graph(6, 6, directed = TRUE)
  iso <- get_isochrone(g, from = "1", lim = c(3, 6, 9))[["1"]]

  expect_true(all(iso[["3"]] %in% iso[["6"]]))
  expect_true(all(iso[["6"]] %in% iso[["9"]]))
  expect_gte(length(iso[["9"]]), length(iso[["6"]]))
  expect_gte(length(iso[["6"]]), length(iso[["3"]]))
})

test_that("setdif = TRUE gives disjoint buckets whose union matches the cumulative bucket", {
  g <- grid_graph(6, 6, directed = TRUE)
  iso_cum <- get_isochrone(g, from = "1", lim = c(3, 6, 9))[["1"]]
  iso_diff <- get_isochrone(g, from = "1", lim = c(3, 6, 9), setdif = TRUE)[["1"]]

  expect_length(intersect(iso_diff[["3"]], iso_diff[["6"]]), 0)
  expect_length(intersect(iso_diff[["6"]], iso_diff[["9"]]), 0)

  union_diff <- Reduce(union, iso_diff)
  expect_setequal(union_diff, iso_cum[["9"]])
})

test_that("the origin node is always included (distance 0)", {
  g <- grid_graph(5, 5, directed = TRUE)
  # a single `lim` value collapses the per-limit nesting level (unlike a
  # vector of limits, see the tests above) - `iso` is directly the node vector.
  iso <- get_isochrone(g, from = "1", lim = 3)[["1"]]
  expect_true("1" %in% iso)
})
