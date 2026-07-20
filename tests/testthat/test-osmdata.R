# Same network shape as test-osmpbf.R's mid-way junction case (node "2" is
# shared between way1's middle and way2's start), built as a fake osmdata
# result instead of a .pbf file: $osm_points carries every node's osm_id
# (including the untagged middle vertex of way1), which is exactly what lets
# graph_from_osmdata() recover the junction that $osm_lines alone could not.
make_test_osmdata <- function(){
  testthat::skip_if_not_installed("sf")

  pts <- sf::st_sf(
    osm_id = c("1", "2", "3", "4"),
    geometry = sf::st_sfc(
      sf::st_point(c(5.00, 47.00)), sf::st_point(c(5.00, 47.01)),
      sf::st_point(c(5.00, 47.02)), sf::st_point(c(5.01, 47.01)),
      crs = 4326)
  )

  lines <- sf::st_sf(
    osm_id = c("101", "102"),
    highway = c("residential", "residential"),
    oneway = c(NA, "yes"),
    geometry = sf::st_sfc(
      sf::st_linestring(matrix(c(5.00,47.00, 5.00,47.01, 5.00,47.02), ncol = 2, byrow = TRUE)),
      sf::st_linestring(matrix(c(5.00,47.01, 5.01,47.01), ncol = 2, byrow = TRUE)),
      crs = 4326)
  )

  structure(list(osm_lines = lines, osm_points = pts), class = "osmdata")
}

test_that("cpp_extract requires the full osmdata object, not just $osm_lines", {
  skip_if_not_installed("sf")
  x <- make_test_osmdata()
  expect_error(cpp_extract(x$osm_lines), "osmdata")
})

test_that("graph_from_osmdata recovers the mid-way junction via $osm_points", {
  skip_if_not_installed("sf")
  x <- make_test_osmdata()

  res <- cppRouting:::graph_from_osmdata(x, output = "data.frame", progress = FALSE)

  expect_equal(nrow(res$coords), 4L)
  # way1 (twoway, split at junction node "2"): A-B + B-C, both directions = 4
  # way2 (oneway=yes): B->D only = 1
  expect_equal(nrow(res$data), 5L)
})

test_that("cpp_extract dispatches a full osmdata object the same way", {
  skip_if_not_installed("sf")
  x <- make_test_osmdata()
  g1 <- cpp_extract(x, output = "data.frame", progress = FALSE)
  g2 <- cppRouting:::graph_from_osmdata(x, output = "data.frame", progress = FALSE)
  expect_equal(g1, g2)
})
