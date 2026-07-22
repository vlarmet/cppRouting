test_shp <- test_path("testdata", "test_network.shp")
graph_from_shapefile <- cppRouting:::graph_from_shapefile

# test_network.prj declares a projected CRS (Lambert-93-like) so that most
# tests below -- which use default coords = TRUE, epsg = NULL -- don't trip
# the geographic-coordinates hard stop. The handful of tests that need a
# genuinely geographic CRS build their own temp copy with this helper,
# instead of mutating the shared fixture.
make_geographic_shp_copy <- function(){
  tmp_base <- tempfile()
  file.copy(test_shp, paste0(tmp_base, ".shp"))
  file.copy(test_path("testdata", "test_network.shx"), paste0(tmp_base, ".shx"))
  file.copy(test_path("testdata", "test_network.dbf"), paste0(tmp_base, ".dbf"))
  writeLines('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]',
             paste0(tmp_base, ".prj"))
  paste0(tmp_base, ".shp")
}

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

test_that("graph_from_shapefile's from/to/node_id are character, not numeric", {
  # A node id is a label, not a quantity -- even Module A's synthetic dense
  # index ("0", "1", ...) is character, matching Module B's real OSM ids and
  # makegraph()'s own dict$ref contract, so cpp_extract() presents one
  # consistent type for node identity regardless of source.
  res <- graph_from_shapefile(test_shp, cost_field = "LONGUEUR", output = "data.frame")
  expect_type(res$data$from, "character")
  expect_type(res$data$to, "character")
  expect_type(res$coords$node_id, "character")
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
  geo_shp <- make_geographic_shp_copy()
  on.exit(unlink(paste0(sub("\\.shp$", "", geo_shp), c(".shp", ".shx", ".dbf", ".prj"))))
  # coords = FALSE isolates this warning (about geometric length reliability)
  # from the separate geographic-coordinates hard stop tested below.
  expect_warning(
    graph_from_shapefile(geo_shp, direction_field = "SENS",
                          direction_values = list(twoway = "Double sens", forward = "Sens unique", backward = "Sens inverse"),
                          snap_digits = 4, output = "data.frame", coords = FALSE),
    "geographic"
  )
})

test_that("graph_from_shapefile stops on geographic coordinates unless coords=FALSE or epsg is set", {
  geo_shp <- make_geographic_shp_copy()
  on.exit(unlink(paste0(sub("\\.shp$", "", geo_shp), c(".shp", ".shx", ".dbf", ".prj"))))
  common <- list(cost_field = "LONGUEUR", direction_field = "SENS",
                  direction_values = list(twoway = "Double sens", forward = "Sens unique", backward = "Sens inverse"),
                  snap_digits = 4, output = "data.frame")

  expect_error(do.call(graph_from_shapefile, c(list(geo_shp), common)), "geographic")

  # coords = FALSE bypasses the issue entirely.
  res <- do.call(graph_from_shapefile, c(list(geo_shp), common, list(coords = FALSE)))
  expect_null(res$coords)

  # epsg reprojects instead of erroring -- the fixture's raw values (0, 100,
  # 200, ..., 500) aren't realistic lon/lat, so just check reprojection
  # actually changed them rather than asserting a specific Lambert-93
  # easting range.
  skip_if_not_installed("sf")
  res2 <- do.call(graph_from_shapefile, c(list(geo_shp), common, list(epsg = 2154)))
  expect_equal(nrow(res2$coords), 7L)
  expect_false(any(res2$coords$X %in% c(0, 100, 200, 300, 400, 500)))
})

test_that("graph_from_shapefile epsg is ignored when the source CRS is already projected", {
  common <- list(cost_field = "LONGUEUR", direction_field = "SENS",
                  direction_values = list(twoway = "Double sens", forward = "Sens unique", backward = "Sens inverse"),
                  snap_digits = 4, output = "data.frame")
  res <- do.call(graph_from_shapefile, c(list(test_shp), common, list(epsg = 2154)))
  res_noepsg <- do.call(graph_from_shapefile, c(list(test_shp), common))
  expect_equal(res$coords, res_noepsg$coords)
})

test_that("graph_from_shapefile output='csv' writes the same aux as output='data.frame' to a separate file", {
  common <- list(cost_field = "LONGUEUR", keep_length = TRUE, keep_fields = "VOCATION", snap_digits = 4)

  out_prefix <- tempfile()
  paths <- do.call(graph_from_shapefile, c(list(test_shp), common, list(output = "csv", out_prefix = out_prefix)))
  on.exit(unlink(paths))
  expect_length(paths, 3L) # edges, coords, aux
  aux_csv <- data.table::fread(paths[3])

  res <- do.call(graph_from_shapefile, c(list(test_shp), common, list(output = "data.frame")))
  expect_equal(as.data.frame(aux_csv), res$aux, ignore_attr = TRUE)
})

test_that("resolve_snap_digits resolves NULL to a CRS-appropriate default, and respects an explicit value", {
  resolve_snap_digits <- cppRouting:::resolve_snap_digits
  expect_equal(resolve_snap_digits(NULL, TRUE), 7L)
  expect_equal(resolve_snap_digits(NULL, FALSE), 4L)
  expect_equal(resolve_snap_digits(3, TRUE), 3)
  expect_equal(resolve_snap_digits(3, FALSE), 3)
})

test_that("graph_from_shapefile's default snap_digits (NULL) picks 7 or 4 depending on the CRS", {
  common <- list(cost_field = "LONGUEUR", direction_field = "SENS",
                  direction_values = list(twoway = "Double sens", forward = "Sens unique", backward = "Sens inverse"),
                  output = "data.frame")

  # Projected fixture: default resolves to 4, same topology as passing it explicitly.
  res_default <- do.call(graph_from_shapefile, c(list(test_shp), common))
  res_explicit4 <- do.call(graph_from_shapefile, c(list(test_shp), common, list(snap_digits = 4)))
  expect_equal(res_default, res_explicit4)

  # Geographic fixture: default resolves to 7, same as passing it explicitly
  # (coords = FALSE dodges the unrelated geographic-coordinates hard stop).
  geo_shp <- make_geographic_shp_copy()
  on.exit(unlink(paste0(sub("\\.shp$", "", geo_shp), c(".shp", ".shx", ".dbf", ".prj"))))
  geo_common <- c(common, list(coords = FALSE))
  res_geo_default <- do.call(graph_from_shapefile, c(list(geo_shp), geo_common))
  res_geo_explicit7 <- do.call(graph_from_shapefile, c(list(geo_shp), geo_common, list(snap_digits = 7)))
  expect_equal(res_geo_default, res_geo_explicit7)
})

test_that("cpp_read_shp_geometry rejects a record whose declared NumPoints exceeds its content length", {
  raw <- readBin(test_shp, "raw", file.info(test_shp)$size)
  # First record starts right after the 100-byte file header; its content
  # (past the 8-byte record header of record-number + content-length, both
  # big-endian int32) starts at byte 108. NumPoints for a PolyLine record
  # sits at content_start + 40, little-endian int32 -- inflate it far beyond
  # what the record's actual content length could hold.
  content_start <- 100L + 8L
  bad <- raw
  bad[(content_start + 40 + 1):(content_start + 40 + 4)] <- as.raw(c(0xFF, 0xFF, 0xFF, 0x7F))
  bad_path <- tempfile(fileext = ".shp")
  on.exit(unlink(bad_path))
  writeBin(bad, bad_path)
  expect_error(cppRouting:::cpp_read_shp_geometry(bad_path), "declares more parts/points")
})

test_that("cpp_read_dbf rejects a record_size too small to hold its own declared fields", {
  dbf_path <- test_path("testdata", "test_network.dbf")
  raw <- readBin(dbf_path, "raw", file.info(dbf_path)$size)
  bad <- raw
  bad[11:12] <- as.raw(c(0x01, 0x00)) # record_size = 1, too small for any real field
  bad_path <- tempfile(fileext = ".dbf")
  on.exit(unlink(bad_path))
  writeBin(bad, bad_path)
  expect_error(cppRouting:::cpp_read_dbf(bad_path), "does not fit within the declared record size")
})

common_shp_args <- list(cost_field = "LONGUEUR", direction_field = "SENS",
                         direction_values = list(twoway = "Double sens", forward = "Sens unique", backward = "Sens inverse"),
                         snap_digits = 4)

test_that("graph_from_shapefile coords=FALSE drops coords from every output mode", {
  res <- do.call(graph_from_shapefile, c(list(test_shp), common_shp_args, list(output = "data.frame", coords = FALSE)))
  expect_null(res$coords)
  expect_equal(nrow(res$data), 10L)

  g <- do.call(graph_from_shapefile, c(list(test_shp), common_shp_args, list(coords = FALSE)))
  expect_null(g$coords)
  expect_equal(g$nbnode, 7L)

  out_prefix <- tempfile()
  paths <- do.call(graph_from_shapefile, c(list(test_shp), common_shp_args,
                                            list(output = "csv", out_prefix = out_prefix, coords = FALSE)))
  on.exit(unlink(paths))
  expect_length(paths, 1L)
  expect_true(file.exists(paths))
  expect_false(file.exists(paste0(out_prefix, "_coords.csv")))
})
