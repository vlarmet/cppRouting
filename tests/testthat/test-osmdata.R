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

  # epsg required: osmdata coordinates are always geographic (WGS84).
  res <- cppRouting:::graph_from_osmdata(x, output = "data.frame", epsg = 2154, progress = FALSE)

  expect_equal(nrow(res$coords), 4L)
  # way1 (twoway, split at junction node "2"): A-B + B-C, both directions = 4
  # way2 (oneway=yes): B->D only = 1
  expect_equal(nrow(res$data), 5L)
})

test_that("cpp_extract dispatches a full osmdata object the same way", {
  skip_if_not_installed("sf")
  x <- make_test_osmdata()
  g1 <- cpp_extract(x, output = "data.frame", coords = FALSE, osm = list(progress = FALSE))
  g2 <- cppRouting:::graph_from_osmdata(x, output = "data.frame", coords = FALSE, progress = FALSE)
  expect_equal(g1, g2)
})

test_that("graph_from_osmdata warns when snap_digits collapses two distinct osm_id points onto the same key", {
  skip_if_not_installed("sf")
  x <- make_test_osmdata()
  # Duplicate node "3"'s coordinate under a different osm_id -- a genuine
  # collision at the default snap_digits precision, which the lookup below
  # can only resolve by silently keeping one of the two ids.
  extra <- sf::st_sf(osm_id = "999",
                      geometry = sf::st_sfc(sf::st_point(c(5.00, 47.02)), crs = 4326))
  x$osm_points <- rbind(x$osm_points, extra)

  # coords = FALSE bypasses the (unrelated) geographic-coordinates hard stop;
  # the collision check itself runs regardless of coords, since it's about
  # topology reconstruction, not the output coords table.
  expect_warning(
    cppRouting:::graph_from_osmdata(x, output = "data.frame", coords = FALSE, progress = FALSE),
    "shared by multiple distinct osm_id"
  )
})

test_that("graph_from_osmdata stops on geographic coordinates unless coords=FALSE or epsg is set", {
  skip_if_not_installed("sf")
  x <- make_test_osmdata()
  expect_error(cppRouting:::graph_from_osmdata(x, output = "data.frame", progress = FALSE), "geographic")
})

test_that("graph_from_osmdata coords=FALSE drops coords from every output mode", {
  skip_if_not_installed("sf")
  x <- make_test_osmdata()

  res <- cppRouting:::graph_from_osmdata(x, output = "data.frame", coords = FALSE, progress = FALSE)
  expect_null(res$coords)
  expect_equal(nrow(res$data), 5L)

  g <- cppRouting:::graph_from_osmdata(x, coords = FALSE, progress = FALSE)
  expect_null(g$coords)

  out_prefix <- tempfile()
  paths <- cppRouting:::graph_from_osmdata(x, output = "csv", out_prefix = out_prefix, coords = FALSE, progress = FALSE)
  on.exit(unlink(paths))
  expect_length(paths, 1L)
  expect_true(file.exists(paths))
  expect_false(file.exists(paste0(out_prefix, "_coords.csv")))
})

test_that("graph_from_osmdata keep_length keeps haversine length (km) as aux, consistent with dist", {
  skip_if_not_installed("sf")
  x <- make_test_osmdata()

  res <- cppRouting:::graph_from_osmdata(x, output = "data.frame", coords = FALSE,
                                          keep_length = TRUE, progress = FALSE)
  expect_false(is.null(res$aux$length))
  # both ways are "residential" (40 km/h, the default profile speed): dist = 1.5*length
  expect_equal(res$data$dist, 1.5 * res$aux$length, tolerance = 1e-6)
})

test_that("graph_from_osmdata output='csv' writes the same aux as output='data.frame' to a separate file", {
  skip_if_not_installed("sf")
  x <- make_test_osmdata()

  out_prefix <- tempfile()
  paths <- cppRouting:::graph_from_osmdata(x, output = "csv", out_prefix = out_prefix, coords = FALSE,
                                            keep_length = TRUE, progress = FALSE)
  on.exit(unlink(paths))
  expect_length(paths, 2L) # edges, aux (no coords)
  aux_csv <- data.table::fread(paths[2])

  res <- cppRouting:::graph_from_osmdata(x, output = "data.frame", coords = FALSE,
                                          keep_length = TRUE, progress = FALSE)
  expect_equal(aux_csv$length, res$aux$length, tolerance = 1e-6)
})
