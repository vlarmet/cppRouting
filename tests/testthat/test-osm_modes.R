# Mode-specific access-tag/oneway resolution (src/io_osm_tags.h resolve_access(),
# resolve_way()) for mode = "foot"/"bicycle"/"hgv" -- added alongside the
# existing "car" (default) coverage in test-osm_tags.R. One tag scenario per
# test, exercised through cpp_extract() on tiny synthetic .pbf fixtures.

make_one_way_pbf <- function(tags){
  path <- tempfile(fileext = ".pbf")
  nodes <- data.frame(id = 1:2, lat = c(47.00, 47.01), lon = c(5.00, 5.00))
  write_test_osmpbf(path, nodes, list(list(refs = c(1, 2), tags = tags)))
  path
}

test_that("profile auto-selection: mode picks the matching default profile when profile is NULL", {
  path <- make_one_way_pbf(list(highway = "footway")) # not in the car/hgv profile, is in foot/bike
  on.exit(unlink(path))

  expect_equal(nrow(cpp_extract(path, output = "data.frame", coords = FALSE, osm = list(progress = FALSE))$data), 0L)
  expect_equal(nrow(cpp_extract(path, output = "data.frame", coords = FALSE, osm = list(mode = "foot", progress = FALSE))$data), 2L)
  expect_equal(nrow(cpp_extract(path, output = "data.frame", coords = FALSE, osm = list(mode = "bicycle", progress = FALSE))$data), 0L) # cycleway, not footway
  expect_equal(nrow(cpp_extract(path, output = "data.frame", coords = FALSE, osm = list(mode = "hgv", progress = FALSE))$data), 0L) # same profile as car
})

test_that("mode = 'foot': foot=yes/foot=no override access, motor_vehicle tags are irrelevant", {
  path1 <- make_one_way_pbf(list(highway = "residential", access = "no", foot = "yes"))
  path2 <- make_one_way_pbf(list(highway = "residential", foot = "no"))
  path3 <- make_one_way_pbf(list(highway = "residential", motor_vehicle = "no")) # foot ignores this
  on.exit(unlink(c(path1, path2, path3)))

  expect_equal(nrow(cpp_extract(path1, output = "data.frame", coords = FALSE, osm = list(mode = "foot", progress = FALSE))$data), 2L)
  expect_equal(nrow(cpp_extract(path2, output = "data.frame", coords = FALSE, osm = list(mode = "foot", progress = FALSE))$data), 0L)
  expect_equal(nrow(cpp_extract(path3, output = "data.frame", coords = FALSE, osm = list(mode = "foot", progress = FALSE))$data), 2L)
})

test_that("mode = 'foot': oneway is always ignored, and use_maxspeed has no effect", {
  path <- make_one_way_pbf(list(highway = "residential", oneway = "yes", maxspeed = "130"))
  on.exit(unlink(path))

  res <- cpp_extract(path, output = "data.frame", coords = FALSE, keep_length = TRUE, osm = list(mode = "foot", progress = FALSE))
  expect_equal(nrow(res$data), 2L) # bidirectional despite oneway=yes

  foot_speed <- osm_foot_profile()$speed[osm_foot_profile()$highway == "residential"]
  expect_equal(res$data$dist, 60 * res$aux$length / foot_speed, tolerance = 1e-6) # maxspeed=130 not applied
})

test_that("mode = 'bicycle': bicycle=yes/bicycle=no override access, vehicle=no also denies", {
  path1 <- make_one_way_pbf(list(highway = "residential", access = "no", bicycle = "yes"))
  path2 <- make_one_way_pbf(list(highway = "residential", bicycle = "no"))
  path3 <- make_one_way_pbf(list(highway = "residential", vehicle = "no"))
  on.exit(unlink(c(path1, path2, path3)))

  expect_equal(nrow(cpp_extract(path1, output = "data.frame", coords = FALSE, osm = list(mode = "bicycle", progress = FALSE))$data), 2L)
  expect_equal(nrow(cpp_extract(path2, output = "data.frame", coords = FALSE, osm = list(mode = "bicycle", progress = FALSE))$data), 0L)
  expect_equal(nrow(cpp_extract(path3, output = "data.frame", coords = FALSE, osm = list(mode = "bicycle", progress = FALSE))$data), 0L)
})

test_that("mode = 'bicycle': oneway:bicycle=no allows contraflow cycling on an otherwise one-way street", {
  path <- make_one_way_pbf(list(highway = "residential", oneway = "yes", `oneway:bicycle` = "no"))
  on.exit(unlink(path))

  res <- cpp_extract(path, output = "data.frame", coords = FALSE, osm = list(mode = "bicycle", progress = FALSE))
  expect_equal(nrow(res$data), 2L) # both directions despite oneway=yes

  # without the override, plain oneway=yes is respected for bicycle same as car
  path2 <- make_one_way_pbf(list(highway = "residential", oneway = "yes"))
  on.exit(unlink(path2), add = TRUE)
  res2 <- cpp_extract(path2, output = "data.frame", coords = FALSE, osm = list(mode = "bicycle", progress = FALSE))
  expect_equal(nrow(res2$data), 1L)
})

test_that("mode = 'bicycle': use_maxspeed has no effect", {
  path <- make_one_way_pbf(list(highway = "residential", maxspeed = "130"))
  on.exit(unlink(path))

  res <- cpp_extract(path, output = "data.frame", coords = FALSE, keep_length = TRUE, osm = list(mode = "bicycle", progress = FALSE))
  bike_speed <- osm_bike_profile()$speed[osm_bike_profile()$highway == "residential"]
  expect_equal(res$data$dist, 60 * res$aux$length / bike_speed, tolerance = 1e-6)
})

test_that("mode = 'hgv': hgv=no drops a way even if otherwise motor-vehicle accessible", {
  path <- make_one_way_pbf(list(highway = "residential", hgv = "no"))
  on.exit(unlink(path))

  expect_equal(nrow(cpp_extract(path, output = "data.frame", coords = FALSE, osm = list(mode = "car", progress = FALSE))$data), 2L)
  expect_equal(nrow(cpp_extract(path, output = "data.frame", coords = FALSE, osm = list(mode = "hgv", progress = FALSE))$data), 0L)
})

test_that("mode = 'hgv': hgv=yes/hgv=destination re-authorize even if access=no", {
  path1 <- make_one_way_pbf(list(highway = "residential", access = "no", hgv = "yes"))
  path2 <- make_one_way_pbf(list(highway = "residential", access = "no", hgv = "destination"))
  on.exit(unlink(c(path1, path2)))

  expect_equal(nrow(cpp_extract(path1, output = "data.frame", coords = FALSE, osm = list(mode = "hgv", progress = FALSE))$data), 2L)
  expect_equal(nrow(cpp_extract(path2, output = "data.frame", coords = FALSE, osm = list(mode = "hgv", progress = FALSE))$data), 2L)
})

test_that("mode = 'hgv': motor_vehicle=no still denies (a general restriction applies to trucks too)", {
  path <- make_one_way_pbf(list(highway = "residential", motor_vehicle = "no"))
  on.exit(unlink(path))

  expect_equal(nrow(cpp_extract(path, output = "data.frame", coords = FALSE, osm = list(mode = "hgv", progress = FALSE))$data), 0L)
})

test_that("mode = 'hgv': use_maxspeed still applies, same as car", {
  path <- make_one_way_pbf(list(highway = "residential", maxspeed = "70"))
  on.exit(unlink(path))

  res <- cpp_extract(path, output = "data.frame", coords = FALSE, keep_length = TRUE, osm = list(mode = "hgv", progress = FALSE))
  expect_equal(res$data$dist, 60 * res$aux$length / 70, tolerance = 1e-6)
})

test_that("osm_foot_profile()/osm_bike_profile() exclude motorway, osm_hgv_profile() equals osm_car_profile()", {
  expect_false("motorway" %in% osm_foot_profile()$highway)
  expect_false("motorway" %in% osm_bike_profile()$highway)
  expect_true("motorway" %in% osm_hgv_profile()$highway)
  expect_equal(osm_hgv_profile(), osm_car_profile())
})

test_that("mode is validated against the known set", {
  path <- make_one_way_pbf(list(highway = "residential"))
  on.exit(unlink(path))
  expect_error(cpp_extract(path, output = "data.frame", osm = list(mode = "unicycle")))
})
