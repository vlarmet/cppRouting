test_shp <- test_path("testdata", "test_network.shp")
graph_from_shapefile <- cppRouting:::graph_from_shapefile
graph_from_sf <- cppRouting:::graph_from_sf

test_that("cpp_extract dispatches .shp paths to the shapefile reader", {
  g1 <- cpp_extract(test_shp, cost_field = "LONGUEUR", direction_field = "SENS",
                     direction_values = list(twoway = "Double sens", forward = "Sens unique", backward = "Sens inverse"),
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
  g1 <- cpp_extract(x, cost_field = "LONGUEUR", direction_field = "SENS",
                     direction_values = list(twoway = "Double sens"))
  g2 <- graph_from_sf(x, cost_field = "LONGUEUR", direction_field = "SENS",
                       direction_values = list(twoway = "Double sens"))
  expect_equal(g1, g2)
})

test_that("cpp_extract rejects unsupported extensions/objects with a clear error", {
  expect_error(cpp_extract("network.gpkg"), "not implemented yet")
  expect_error(cpp_extract("network.pbf"), "File not found") # .pbf IS implemented, this one just doesn't exist
  expect_error(cpp_extract("network.osm"), "Plain OSM XML")
  expect_error(cpp_extract("network.txt"), "Unrecognized file extension")
  expect_error(cpp_extract(list(a = 1)), "does not know how to handle")
})

test_that("cpp_extract refuses OSM-derived objects instead of silently mishandling them", {
  fake_osmdata <- structure(list(osm_lines = NULL), class = "osmdata")
  expect_error(cpp_extract(fake_osmdata), "OSM-derived")

  skip_if_not_installed("sf")
  x <- sf::st_sf(osm_id = "123", highway = "residential",
                 geometry = sf::st_sfc(sf::st_linestring(matrix(c(0,0,100,0), ncol = 2, byrow = TRUE)), crs = 2154))
  expect_error(cpp_extract(x), "OSM-derived")
})
