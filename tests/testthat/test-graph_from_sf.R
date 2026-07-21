# Same 5-feature fixture as test-graph_from_shapefile.R (see there for the
# hand-derived expected topology: 7 nodes, 10 directed edges), built directly
# as an sf object instead of read from a file, to exercise graph_from_sf()'s
# own st_coordinates()-based extraction path.

graph_from_sf <- cppRouting:::graph_from_sf

make_test_sf <- function(){
  testthat::skip_if_not_installed("sf")
  geoms <- sf::st_sfc(
    sf::st_linestring(matrix(c(0,0, 100,0), ncol = 2, byrow = TRUE)),
    sf::st_linestring(matrix(c(100,0, 200,0), ncol = 2, byrow = TRUE)),
    sf::st_linestring(matrix(c(200.00000001,0.00000001, 300,0), ncol = 2, byrow = TRUE)),
    sf::st_multilinestring(list(
      matrix(c(300,0, 400,0), ncol = 2, byrow = TRUE),
      matrix(c(400,0, 400,100), ncol = 2, byrow = TRUE)
    )),
    sf::st_linestring(matrix(c(400.00000002,100.00000003, 500,100), ncol = 2, byrow = TRUE)),
    crs = 2154 # RGF93 / Lambert-93, projected (metres) -- same family as ROUTE500
  )
  sf::st_sf(
    VOCATION = c("Liaison locale", "Liaison principale", "Type autoroutier", "Liaison locale", "Liaison locale"),
    ACCES = c(NA, NA, "Libre", NA, NA),
    SENS = c("Double sens", "Sens unique", "Sens inverse", "Double sens", "Double sens"),
    LONGUEUR = c(100, 100, 100, 141.42, 100),
    geometry = geoms
  )
}

test_that("graph_from_sf builds the expected topology and direction expansion", {
  x <- make_test_sf()
  res <- graph_from_sf(x,
                        cost_field = "LONGUEUR",
                        direction_field = "SENS",
                        direction_values = list(twoway = "Double sens", forward = "Sens unique", backward = "Sens inverse"),
                        snap_digits = 4,
                        output = "data.frame")

  expect_equal(nrow(res$coords), 7L)
  expect_equal(nrow(res$data), 10L)

  key <- function(x, y) sprintf("%.2f_%.2f", round(x, 4), round(y, 4))
  id2key <- setNames(key(res$coords$X, res$coords$Y), res$coords$node_id)
  got <- data.frame(from_key = id2key[as.character(res$data$from)],
                     to_key = id2key[as.character(res$data$to)],
                     dist = res$data$dist, row.names = NULL)

  expected <- data.frame(
    from_key = c(key(0,0), key(100,0), key(100,0), key(300,0), key(300,0),
                 key(400,0), key(400,0), key(400,100), key(400,100), key(500,100)),
    to_key   = c(key(100,0), key(0,0), key(200,0), key(200,0), key(400,0),
                 key(300,0), key(400,100), key(400,0), key(500,100), key(400,100)),
    dist     = c(100, 100, 100, 100, 141.42, 141.42, 141.42, 141.42, 100, 100)
  )

  merged <- merge(got, expected, by = c("from_key", "to_key"), suffixes = c("_got", "_exp"))
  expect_equal(nrow(merged), 10L)
  expect_equal(merged$dist_got, merged$dist_exp, tolerance = 1e-6)
})

test_that("graph_from_sf output='graph' returns a usable graph", {
  x <- make_test_sf()
  g <- graph_from_sf(x,
                      cost_field = "LONGUEUR",
                      direction_field = "SENS",
                      direction_values = list(twoway = "Double sens", forward = "Sens unique", backward = "Sens inverse"),
                      snap_digits = 4)
  expect_equal(g$nbnode, 7L)
  d <- get_distance_pair(g, from = g$dict$ref[1], to = g$dict$ref[nrow(g$dict)])
  expect_equal(d, 200)
})

test_that("graph_from_sf speed_field without length_field converts geometric length (meters) to km", {
  # Regression test: geom_length is computed directly from the (projected,
  # meter-scale) coordinates -- without converting to km first, 60*length/speed
  # silently produced a cost 1000x too large (e.g. a 100m segment at 60km/h
  # should cost 0.1 minutes, not 100).
  x <- make_test_sf()
  x$VITESSE <- 60
  res <- graph_from_sf(x, speed_field = "VITESSE", directed = TRUE, snap_digits = 4, output = "data.frame")
  # First segment: (0,0) -> (100,0), 100m = 0.1km, at 60km/h -> 60*0.1/60 = 0.1 min
  expect_equal(min(res$data$dist), 0.1, tolerance = 1e-6)
})

test_that("graph_from_sf speed_field without length_field does NOT convert for geographic coordinates", {
  testthat::skip_if_not_installed("sf")
  x <- sf::st_transform(make_test_sf(), 4326)
  x$VITESSE <- 60
  res <- graph_from_sf(x, speed_field = "VITESSE", directed = TRUE, snap_digits = 4,
                        output = "data.frame", coords = FALSE)
  # Geometric length here is in degrees (unreliable regardless, already
  # warned about by warn_if_geographic_without_length) -- just confirm no
  # /1000 conversion silently applies to a non-meter unit.
  expect_true(all(res$data$dist > 0))
})

test_that("graph_from_sf keep_length keeps geometric length (km) as aux, independent of cost_field", {
  testthat::skip_if_not_installed("sf")
  x <- sf::st_sf(
    id = 1:2,
    cost = c(999, 999), # arbitrary, unrelated to geometry -- proves keep_length doesn't depend on it
    geometry = sf::st_sfc(
      sf::st_linestring(matrix(c(0,0, 100,0), ncol = 2, byrow = TRUE)),     # 100m = 0.1km
      sf::st_linestring(matrix(c(100,0, 100,300), ncol = 2, byrow = TRUE)), # 300m = 0.3km
      crs = 2154
    )
  )
  res <- graph_from_sf(x, cost_field = "cost", keep_length = TRUE, directed = TRUE, output = "data.frame")
  expect_equal(sort(res$aux$length), c(0.1, 0.3))
  expect_true(all(res$data$dist == 999)) # cost_field unaffected by keep_length

  # combined with keep_fields: both aux columns present
  res2 <- graph_from_sf(x, cost_field = "cost", keep_fields = "id", keep_length = TRUE,
                         directed = TRUE, output = "data.frame")
  expect_setequal(names(res2$aux), c("id", "length"))
})

test_that("graph_from_sf output='csv' writes the same aux as output='data.frame' to a separate file", {
  testthat::skip_if_not_installed("sf")
  x <- make_test_sf()
  common <- list(cost_field = "LONGUEUR", keep_length = TRUE, snap_digits = 4)

  out_prefix <- tempfile()
  paths <- do.call(graph_from_sf, c(list(x), common, list(output = "csv", out_prefix = out_prefix)))
  on.exit(unlink(paths))
  expect_length(paths, 3L) # edges, coords, aux
  aux_csv <- data.table::fread(paths[3])

  res <- do.call(graph_from_sf, c(list(x), common, list(output = "data.frame")))
  expect_equal(aux_csv$length, res$aux$length)
})

test_that("graph_from_sf rejects non-sf input and non-linestring geometries", {
  testthat::skip_if_not_installed("sf")
  expect_error(graph_from_sf(data.frame(x = 1)), "sf object")

  pts <- sf::st_sf(id = 1, geometry = sf::st_sfc(sf::st_point(c(0, 0)), crs = 2154))
  expect_error(graph_from_sf(pts), "LINESTRING")
})

test_that("graph_from_sf warns on geographic CRS without an explicit length/cost field", {
  testthat::skip_if_not_installed("sf")
  x <- make_test_sf()
  x <- sf::st_transform(x, 4326)
  # coords = FALSE isolates this warning (about geometric length reliability)
  # from the separate geographic-coordinates hard stop tested below.
  expect_warning(
    graph_from_sf(x, direction_field = "SENS",
                  direction_values = list(twoway = "Double sens", forward = "Sens unique", backward = "Sens inverse"),
                  snap_digits = 4, output = "data.frame", coords = FALSE),
    "geographic"
  )
})

test_that("graph_from_sf stops on geographic coordinates unless coords=FALSE or epsg is set", {
  testthat::skip_if_not_installed("sf")
  x <- sf::st_transform(make_test_sf(), 4326)
  common <- list(cost_field = "LONGUEUR", direction_field = "SENS",
                  direction_values = list(twoway = "Double sens", forward = "Sens unique", backward = "Sens inverse"),
                  snap_digits = 4, output = "data.frame")

  expect_error(do.call(graph_from_sf, c(list(x), common)), "geographic")

  res <- do.call(graph_from_sf, c(list(x), common, list(coords = FALSE)))
  expect_null(res$coords)

  # Round-tripping through 4326 and back to the original projected CRS
  # (2154) should approximately recover the original small-number Lambert-93
  # coordinates (0..500) used to build the fixture -- not a coincidence with
  # a large real-world easting, since this fixture isn't placed in France.
  res2 <- do.call(graph_from_sf, c(list(x), common, list(epsg = 2154)))
  expect_equal(nrow(res2$coords), 7L)
  expect_true(all(res2$coords$X > -10 & res2$coords$X < 510))
})

test_that("graph_from_sf epsg is ignored when the source CRS is already projected", {
  x <- make_test_sf() # crs = 2154, already projected
  common <- list(cost_field = "LONGUEUR", direction_field = "SENS",
                  direction_values = list(twoway = "Double sens", forward = "Sens unique", backward = "Sens inverse"),
                  snap_digits = 4, output = "data.frame")
  res <- do.call(graph_from_sf, c(list(x), common, list(epsg = 3857)))
  res_noepsg <- do.call(graph_from_sf, c(list(x), common))
  expect_equal(res$coords, res_noepsg$coords)
})
