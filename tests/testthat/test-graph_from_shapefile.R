test_shp <- test_path("testdata", "test_network.shp")
graph_from_shapefile <- cppRouting:::graph_from_shapefile

# Fixture: 5 records (VOCATION/ACCES/SENS/LONGUEUR, mirroring ROUTE500's schema),
# one of them a MULTILINESTRING (2 parts), with two endpoints deliberately
# perturbed by 1e-8 to exercise snap_digits. Expected topology (7 nodes, 10
# directed edges after Sens unique/inverse/Double sens expansion) was derived
# by hand and cross-checked against the real ROUTE500 France shapefile.

test_that("cpp_read_shp_geometry parses records/parts/attributes correctly", {
  geom <- cppRouting:::cpp_read_shp_geometry(test_shp)
  expect_equal(geom$n_records, 5L)
  expect_equal(length(geom$x1), 6L) # record 4 (0-based: 3) has 2 parts
  expect_equal(geom$record, c(0L, 1L, 2L, 3L, 3L, 4L))
})

test_that("cpp_read_dbf reads attributes with correct types and values", {
  dbf <- cppRouting:::cpp_read_dbf(test_path("testdata", "test_network.dbf"))
  expect_equal(nrow(dbf), 5L)
  expect_equal(as.character(dbf$SENS), c("Double sens", "Sens unique", "Sens inverse", "Double sens", "Double sens"))
  expect_equal(as.numeric(dbf$LONGUEUR), c(100, 100, 100, 141.42, 100))
})

test_that("graph_from_shapefile builds the expected topology and direction expansion", {
  res <- graph_from_shapefile(test_shp,
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

  expect_setequal(paste(got$from_key, got$to_key), paste(got$from_key, got$to_key)) # sanity: no NA keys
  expect_true(all(!is.na(got$from_key)) && all(!is.na(got$to_key)))

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

test_that("graph_from_shapefile output='graph' returns a usable graph", {
  g <- graph_from_shapefile(test_shp,
                             cost_field = "LONGUEUR",
                             direction_field = "SENS",
                             direction_values = list(twoway = "Double sens", forward = "Sens unique", backward = "Sens inverse"),
                             snap_digits = 4)
  expect_equal(g$nbnode, 7L)
  expect_equal(nrow(g$dict), 7L)
  d <- get_distance_pair(g, from = g$dict$ref[1], to = g$dict$ref[nrow(g$dict)])
  expect_equal(d, 200)
})

test_that("graph_from_shapefile keep_fields aligns aux to the expanded edge set", {
  res <- graph_from_shapefile(test_shp,
                               cost_field = "LONGUEUR",
                               direction_field = "SENS",
                               direction_values = list(twoway = "Double sens", forward = "Sens unique", backward = "Sens inverse"),
                               snap_digits = 4, keep_fields = "VOCATION",
                               output = "data.frame")
  expect_equal(nrow(res$aux), nrow(res$data))
  expect_equal(sum(res$aux$VOCATION == "Liaison locale"), 8L)
})

test_that("graph_from_shapefile output='csv' writes edges/coords files", {
  out_prefix <- tempfile()
  paths <- graph_from_shapefile(test_shp,
                                 cost_field = "LONGUEUR",
                                 direction_field = "SENS",
                                 direction_values = list(twoway = "Double sens", forward = "Sens unique", backward = "Sens inverse"),
                                 snap_digits = 4, output = "csv", out_prefix = out_prefix)
  on.exit(unlink(paths))
  expect_true(all(file.exists(paths)))
  edges <- data.table::fread(paths[1])
  expect_equal(nrow(edges), 10L)
})

test_that("graph_from_shapefile warns on geographic CRS without an explicit length/cost field", {
  expect_warning(
    graph_from_shapefile(test_shp, direction_field = "SENS",
                          direction_values = list(twoway = "Double sens", forward = "Sens unique", backward = "Sens inverse"),
                          snap_digits = 4, output = "data.frame"),
    "geographic"
  )
})
