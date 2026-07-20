test_that("cpp_osmpbf_counts reads back the exact synthetic fixture", {
  path <- make_test_pbf()
  on.exit(unlink(path))
  counts <- cppRouting:::cpp_osmpbf_counts(path)
  expect_equal(counts$n_nodes, 4)
  expect_equal(counts$n_ways, 2)
  expect_equal(counts$n_tagged_ways, 2)
})

test_that("cpp_extract on a .pbf splits a way at a mid-way junction (node 2)", {
  path <- make_test_pbf()
  on.exit(unlink(path))

  res <- cpp_extract(path, output = "data.frame", progress = FALSE)

  # 4 osm nodes in, but only 3 become graph nodes: node 2 (the junction),
  # and the two other endpoints referenced only once each collapse per way
  # end -- node 1, node 2, node 3, node 4 are ALL split points here (1 and 3
  # are way1's own endpoints, 4 is way2's endpoint, 2 is shared) so all 4
  # remain as graph nodes.
  expect_equal(nrow(res$coords), 4L)

  # way A-B-C (residential, twoway by default, no oneway tag) split at B (node 2):
  #   edge A-B (both directions) + edge B-C (both directions) = 4 directed edges
  # way B-D (oneway=yes): B->D only = 1 directed edge
  expect_equal(nrow(res$data), 5L)
})

test_that("cpp_extract on a .pbf resolves oneway correctly and drops non-profile highways", {
  path <- tempfile(fileext = ".pbf")
  on.exit(unlink(path))
  nodes <- data.frame(id = 1:2, lat = c(47.00, 47.01), lon = c(5.00, 5.00))
  ways <- list(list(refs = c(1, 2), tags = list(highway = "footway"))) # not in the default car profile
  write_test_osmpbf(path, nodes, ways)

  res <- cpp_extract(path, output = "data.frame", progress = FALSE)
  expect_equal(nrow(res$data), 0L)
  expect_equal(nrow(res$coords), 0L)
})
