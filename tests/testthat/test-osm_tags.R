# Tag-resolution logic (src/io_osm_tags.h resolve_way()): highway->profile
# lookup, the access/motor_vehicle/motorcar/vehicle allow-deny cascade,
# maxspeed override (use_maxspeed), and oneway resolution (explicit tag,
# assume_oneway_motorway inference). Exercised through cpp_extract() on tiny
# synthetic .pbf fixtures built with write_test_osmpbf() (see
# helper-osmpbf.R), one tag scenario per test -- none of this had any
# regression coverage before, despite being the core business logic that
# decides cost/direction for every OSM-sourced edge.

profile <- osm_car_profile()
residential_speed <- profile$speed[profile$highway == "residential"]

make_one_way_pbf <- function(tags){
  path <- tempfile(fileext = ".pbf")
  nodes <- data.frame(id = 1:2, lat = c(47.00, 47.01), lon = c(5.00, 5.00))
  write_test_osmpbf(path, nodes, list(list(refs = c(1, 2), tags = tags)))
  path
}

test_that("use_maxspeed overrides the profile default speed (plain km/h)", {
  path <- make_one_way_pbf(list(highway = "residential", maxspeed = "70"))
  on.exit(unlink(path))

  res <- cpp_extract(path, output = "data.frame", coords = FALSE, keep_length = TRUE, osm = list(progress = FALSE))
  expect_equal(res$data$dist, 60 * res$aux$length / 70, tolerance = 1e-6)

  res_off <- cpp_extract(path, output = "data.frame", coords = FALSE, keep_length = TRUE,
                          osm = list(use_maxspeed = FALSE, progress = FALSE))
  expect_equal(res_off$data$dist, 60 * res_off$aux$length / residential_speed, tolerance = 1e-6)
})

test_that("use_maxspeed parses a '<n> mph' value", {
  path <- make_one_way_pbf(list(highway = "residential", maxspeed = "30 mph"))
  on.exit(unlink(path))

  res <- cpp_extract(path, output = "data.frame", coords = FALSE, keep_length = TRUE, osm = list(progress = FALSE))
  expect_equal(res$data$dist, 60 * res$aux$length / (30 * 1.609344), tolerance = 1e-6)
})

test_that("an unparseable maxspeed value falls back to the profile default speed", {
  path <- make_one_way_pbf(list(highway = "residential", maxspeed = "FR:urban"))
  on.exit(unlink(path))

  res <- cpp_extract(path, output = "data.frame", coords = FALSE, keep_length = TRUE, osm = list(progress = FALSE))
  expect_equal(res$data$dist, 60 * res$aux$length / residential_speed, tolerance = 1e-6)
})

test_that("assume_oneway_motorway makes an untagged motorway one-way by default", {
  path <- make_one_way_pbf(list(highway = "motorway")) # no oneway tag
  on.exit(unlink(path))

  res <- cpp_extract(path, output = "data.frame", coords = FALSE, osm = list(progress = FALSE))
  expect_equal(nrow(res$data), 1L) # forward only

  res_off <- cpp_extract(path, output = "data.frame", coords = FALSE,
                          osm = list(assume_oneway_motorway = FALSE, progress = FALSE))
  expect_equal(nrow(res_off$data), 2L) # twoway
})

test_that("an explicit oneway=no on a motorway overrides assume_oneway_motorway", {
  path <- make_one_way_pbf(list(highway = "motorway", oneway = "no"))
  on.exit(unlink(path))

  res <- cpp_extract(path, output = "data.frame", coords = FALSE, osm = list(progress = FALSE)) # assume_oneway_motorway=TRUE (default)
  expect_equal(nrow(res$data), 2L) # explicit tag wins over the motorway inference
})

test_that("oneway=-1 and oneway=reverse both produce a single backward-only edge", {
  path1 <- make_one_way_pbf(list(highway = "residential", oneway = "-1"))
  path2 <- make_one_way_pbf(list(highway = "residential", oneway = "reverse"))
  on.exit(unlink(c(path1, path2)))

  res1 <- cpp_extract(path1, output = "data.frame", coords = FALSE, osm = list(progress = FALSE))
  res2 <- cpp_extract(path2, output = "data.frame", coords = FALSE, osm = list(progress = FALSE))
  expect_equal(nrow(res1$data), 1L)
  expect_equal(res1$data[, c("from", "to")], res2$data[, c("from", "to")])
  # backward: from the way's last node (osm id 2) to its first node (osm id
  # 1) -- from/to carry the real OSM node id (as character), not a remapped
  # dense index.
  expect_equal(res1$data$from, "2")
  expect_equal(res1$data$to, "1")
})

test_that("access=no or access=private drops a way", {
  path1 <- make_one_way_pbf(list(highway = "residential", access = "no"))
  path2 <- make_one_way_pbf(list(highway = "residential", access = "private"))
  on.exit(unlink(c(path1, path2)))

  expect_equal(nrow(cpp_extract(path1, output = "data.frame", coords = FALSE, osm = list(progress = FALSE))$data), 0L)
  expect_equal(nrow(cpp_extract(path2, output = "data.frame", coords = FALSE, osm = list(progress = FALSE))$data), 0L)
})

test_that("motor_vehicle=yes / motorcar=yes override access=no", {
  path1 <- make_one_way_pbf(list(highway = "residential", access = "no", motor_vehicle = "yes"))
  path2 <- make_one_way_pbf(list(highway = "residential", access = "no", motorcar = "yes"))
  on.exit(unlink(c(path1, path2)))

  # no oneway tag -> twoway by default (2 directed edges), kept rather than dropped
  expect_equal(nrow(cpp_extract(path1, output = "data.frame", coords = FALSE, osm = list(progress = FALSE))$data), 2L)
  expect_equal(nrow(cpp_extract(path2, output = "data.frame", coords = FALSE, osm = list(progress = FALSE))$data), 2L)
})

test_that("motor_vehicle=no / vehicle=no drop a way even without access=no", {
  path1 <- make_one_way_pbf(list(highway = "residential", motor_vehicle = "no"))
  path2 <- make_one_way_pbf(list(highway = "residential", vehicle = "no"))
  on.exit(unlink(c(path1, path2)))

  expect_equal(nrow(cpp_extract(path1, output = "data.frame", coords = FALSE, osm = list(progress = FALSE))$data), 0L)
  expect_equal(nrow(cpp_extract(path2, output = "data.frame", coords = FALSE, osm = list(progress = FALSE))$data), 0L)
})

test_that("a custom profile table is honored instead of the default osm_car_profile()", {
  path <- make_one_way_pbf(list(highway = "track")) # not in the default car profile
  on.exit(unlink(path))

  res_default <- cpp_extract(path, output = "data.frame", coords = FALSE, osm = list(progress = FALSE))
  expect_equal(nrow(res_default$data), 0L)

  custom_profile <- data.frame(highway = "track", clazz = 99L, speed = 15)
  res_custom <- cpp_extract(path, output = "data.frame", coords = FALSE, keep_length = TRUE,
                             osm = list(profile = custom_profile, progress = FALSE))
  expect_equal(nrow(res_custom$data), 2L) # no oneway tag -> twoway by default
  expect_equal(res_custom$data$dist, 60 * res_custom$aux$length / 15, tolerance = 1e-6)
})
