test_shp <- test_path("testdata", "test_network.shp")
graph_from_shapefile <- cppRouting:::graph_from_shapefile
graph_from_sf <- cppRouting:::graph_from_sf

test_that("cpp_extract dispatches .shp paths to the shapefile reader", {
  g1 <- cpp_extract(test_shp,
                     geom_params = list(cost_field = "LONGUEUR", direction_field = "SENS",
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
  g1 <- cpp_extract(x, geom_params = list(cost_field = "LONGUEUR", direction_field = "SENS",
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

test_that("cpp_extract's osm_params/geom_params list arguments reject unrecognized element names", {
  # A flat argument typo (e.g. cpp_extract(x, us_maxspeed = FALSE)) is caught
  # for free by R's own "unused argument" error; a typo inside a list
  # (osm_params = list(us_maxspeed = FALSE)) would otherwise silently keep the
  # default instead -- merge_arg_list() must catch this explicitly.
  expect_error(cpp_extract(test_shp, geom_params = list(cost_feild = "LONGUEUR")), "unrecognized element")
  expect_error(cpp_extract("dummy.pbf", osm_params = list(us_maxspeed = FALSE)), "unrecognized element")
})

test_that("cpp_extract's osm_params/geom_params list arguments must actually be lists", {
  expect_error(cpp_extract(test_shp, geom_params = "LONGUEUR"), "must be a list")
  expect_error(cpp_extract("dummy.pbf", osm_params = TRUE), "must be a list")
})

test_that("cpp_extract's osm_params/geom_params accept a partial list, keeping unspecified defaults", {
  g1 <- cpp_extract(test_shp, geom_params = list(cost_field = "LONGUEUR"))
  g2 <- cpp_extract(test_shp,
                     geom_params = list(cost_field = "LONGUEUR", speed_field = NULL, length_field = NULL,
                                     direction_field = NULL, direction_values = NULL,
                                     directed = FALSE, keep_fields = NULL))
  expect_equal(g1, g2)
})

test_that("output = 'graph' carries a single-column aux (keep_fields) as a named data.frame", {
  # makegraph()'s aux natively accepts a named data.frame (one or more
  # numeric columns, each an additional aggregate-summable weight -- see
  # get_distance_pair(..., aggregate_aux = TRUE)); cpp_extract()'s aux is
  # passed straight through, column names preserved.
  g <- cpp_extract(test_shp, geom_params = list(cost_field = "LONGUEUR", keep_fields = "LONGUEUR"), output = "graph")
  expect_true(is.data.frame(g$attrib$aux))
  expect_equal(names(g$attrib$aux), "LONGUEUR")
  expect_equal(g$attrib$aux$LONGUEUR, g$data$dist) # LONGUEUR was also used as cost_field here
})

test_that("output = 'graph' carries a multi-column numeric aux (keep_fields + keep_length) intact", {
  g <- cpp_extract(test_shp, geom_params = list(cost_field = "LONGUEUR", keep_fields = "LONGUEUR"),
                    keep_length = TRUE, output = "graph")
  expect_true(is.data.frame(g$attrib$aux))
  expect_setequal(names(g$attrib$aux), c("LONGUEUR", "length"))
  expect_equal(g$attrib$aux$LONGUEUR, g$data$dist)
})

test_that("output = 'graph' with a non-numeric keep_fields column errors via makegraph()'s aux validation", {
  # keep_fields accepts any column type for output = "data.frame"/"csv" (it's
  # purely informational there), but output = "graph" routes it through
  # makegraph()'s aux, which requires every column to be numeric (it's an
  # aggregate-summable edge weight) -- a character column like SENS
  # correctly errors instead of silently doing something meaningless with it.
  expect_error(
    cpp_extract(test_shp, geom_params = list(cost_field = "LONGUEUR", keep_fields = c("LONGUEUR", "SENS")), output = "graph"),
    "must be numeric"
  )
})
