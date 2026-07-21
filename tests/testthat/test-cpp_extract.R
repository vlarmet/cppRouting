test_shp <- test_path("testdata", "test_network.shp")
graph_from_shapefile <- cppRouting:::graph_from_shapefile
graph_from_sf <- cppRouting:::graph_from_sf

test_that("cpp_extract dispatches .shp paths to the shapefile reader", {
  g1 <- cpp_extract(test_shp,
                     non_osm = list(cost_field = "LONGUEUR", direction_field = "SENS",
                                     direction_values = list(twoway = "Double sens", forward = "Sens unique", backward = "Sens inverse")),
                     snap_digits = 4)
  g2 <- graph_from_shapefile(test_shp, cost_field = "LONGUEUR", direction_field = "SENS",
                              direction_values = list(twoway = "Double sens", forward = "Sens unique", backward = "Sens inverse"),
                              snap_digits = 4)
  expect_equal(g1, g2)
})

test_that("cpp_extract dispatches sf objects to the sf reader", {
  skip_if_not_installed("sf")
  x <- sf::st_sf(SENS = c("Double sens"), LONGUEUR = 100,
                 geometry = sf::st_sfc(sf::st_linestring(matrix(c(0,0,100,0), ncol = 2, byrow = TRUE)), crs = 2154))
  g1 <- cpp_extract(x, non_osm = list(cost_field = "LONGUEUR", direction_field = "SENS",
                                        direction_values = list(twoway = "Double sens")))
  g2 <- graph_from_sf(x, cost_field = "LONGUEUR", direction_field = "SENS",
                       direction_values = list(twoway = "Double sens"))
  expect_equal(g1, g2)
})

test_that("cpp_extract rejects unsupported extensions/objects with a clear error", {
  expect_error(cpp_extract("network.gpkg"), "not supported")
  # coords = FALSE: the geographic-coordinates hard stop is checked before
  # file existence, and would otherwise mask the "File not found" this test
  # targets (.pbf IS implemented, this path just doesn't exist).
  expect_error(cpp_extract("network.pbf", coords = FALSE), "File not found")
  expect_error(cpp_extract("network.osm"), "Plain OSM XML")
  expect_error(cpp_extract("network.txt"), "Unrecognized file extension")
  expect_error(cpp_extract(list(a = 1)), "does not know how to handle")
})

test_that("cpp_extract refuses OSM-derived objects it can't safely handle", {
  fake_osmdata <- structure(list(osm_lines = NULL), class = "osmdata")
  expect_error(cpp_extract(fake_osmdata), "osm_points")

  skip_if_not_installed("sf")
  x <- sf::st_sf(osm_id = "123", highway = "residential",
                 geometry = sf::st_sfc(sf::st_linestring(matrix(c(0,0,100,0), ncol = 2, byrow = TRUE)), crs = 2154))
  expect_error(cpp_extract(x), "osmdata")
})

test_that("cpp_extract's osm/non_osm list arguments reject unrecognized element names", {
  # A flat argument typo (e.g. cpp_extract(x, us_maxspeed = FALSE)) is caught
  # for free by R's own "unused argument" error; a typo inside a list
  # (osm = list(us_maxspeed = FALSE)) would otherwise silently keep the
  # default instead -- merge_arg_list() must catch this explicitly.
  expect_error(cpp_extract(test_shp, non_osm = list(cost_feild = "LONGUEUR")), "unrecognized element")
  expect_error(cpp_extract("dummy.pbf", osm = list(us_maxspeed = FALSE)), "unrecognized element")
})

test_that("cpp_extract's osm/non_osm list arguments must actually be lists", {
  expect_error(cpp_extract(test_shp, non_osm = "LONGUEUR"), "must be a list")
  expect_error(cpp_extract("dummy.pbf", osm = TRUE), "must be a list")
})

test_that("cpp_extract's osm/non_osm accept a partial list, keeping unspecified defaults", {
  g1 <- cpp_extract(test_shp, non_osm = list(cost_field = "LONGUEUR"))
  g2 <- cpp_extract(test_shp,
                     non_osm = list(cost_field = "LONGUEUR", speed_field = NULL, length_field = NULL,
                                     direction_field = NULL, direction_values = NULL,
                                     directed = FALSE, keep_fields = NULL))
  expect_equal(g1, g2)
})
