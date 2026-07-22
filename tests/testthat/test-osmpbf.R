test_that("cpp_osmpbf_counts reads back the exact synthetic fixture", {
  path <- make_test_pbf()
  on.exit(unlink(path))
  counts <- cppRouting:::cpp_osmpbf_counts(path)
  expect_equal(counts$n_nodes, 4)
  expect_equal(counts$n_ways, 2)
  expect_equal(counts$n_tagged_ways, 2)
})

test_that("cpp_extract on a .pbf splits a way at a mid-way junction (node 2)", {
  skip_if_not_installed("sf")
  path <- make_test_pbf()
  on.exit(unlink(path))

  # epsg required: OSM coordinates are always geographic (WGS84).
  res <- cpp_extract(path, output = "data.frame", epsg = 2154, osm = list(progress = FALSE))

  # 4 osm nodes in, but only 3 become graph nodes: node 2 (the junction),
  # and the two other endpoints referenced only once each collapse per way
  # end -- node 1, node 2, node 3, node 4 are ALL split points here (1 and 3
  # are way1's own endpoints, 4 is way2's endpoint, 2 is shared) so all 4
  # remain as graph nodes.
  expect_equal(nrow(res$coords), 4L)

  # way A-B-C (residential, twoway by default, no oneway tag) split at B (node 2):
  #   edge A-B (both directions) + edge B-C (both directions) = 4 directed edges
  # way B-D (oneway=yes): B->D only = 1 directed edge
  expect_equal(nrow(res$data), 5L)
})

test_that("cpp_extract on a .pbf keeps the real OSM node id, not a remapped dense index", {
  path <- make_test_pbf()
  on.exit(unlink(path))

  res <- cpp_extract(path, output = "data.frame", coords = FALSE, osm = list(progress = FALSE))
  # nodes 1-4 in the fixture (see helper-osmpbf.R); from/to must be a subset
  # of the real ids, not a synthetic 0-based renumbering. Character, not
  # numeric: a node id is a label, not a quantity (matches makegraph()'s own
  # dict$ref contract).
  expect_type(res$data$from, "character")
  expect_setequal(c(res$data$from, res$data$to), 1:4)

  full <- cpp_extract(path, output = "data.frame", epsg = 2154, osm = list(progress = FALSE))
  expect_type(full$coords$node_id, "character")
  expect_setequal(full$coords$node_id, 1:4)
})

test_that("a large OSM id survives epsg reprojection's csv rewrite pass without integer64 corruption", {
  skip_if_not_installed("sf")
  path <- tempfile(fileext = ".pbf")
  on.exit(unlink(path))
  big_id <- 6564874389 # a real-looking OSM node id, well past 2^31 (~2.147e9)
  nodes <- data.frame(id = c(1, big_id), lat = c(47.00, 47.01), lon = c(5.00, 5.00))
  ways <- list(list(refs = c(1, big_id), tags = list(highway = "residential")))
  write_test_osmpbf(path, nodes, ways)

  # epsg triggers reproject_coords_csv(): it reads coords.csv back with
  # fread() (safe, since the C++ writer already appended ".0" to big_id),
  # reprojects X/Y, then rewrites the file with fwrite() -- which, left
  # unguarded, would print the now-in-memory whole-number double back out as
  # plain digits, silently re-introducing the exact ambiguity the ".0" was
  # there to prevent.
  out_prefix <- tempfile()
  paths <- cpp_extract(path, output = "csv", out_prefix = out_prefix, epsg = 2154, osm = list(progress = FALSE))
  on.exit(unlink(paths), add = TRUE)
  coords_csv <- data.table::fread(paths[2])
  expect_type(coords_csv$node_id, "double") # not bit64::integer64
  expect_true(big_id %in% coords_csv$node_id)
})

test_that("a real OSM id past the 32-bit integer range round-trips through csv without integer64 corruption", {
  path <- tempfile(fileext = ".pbf")
  on.exit(unlink(path))
  big_id <- 6564874389 # a real-looking OSM node id, well past 2^31 (~2.147e9)
  nodes <- data.frame(id = c(1, big_id), lat = c(47.00, 47.01), lon = c(5.00, 5.00))
  ways <- list(list(refs = c(1, big_id), tags = list(highway = "residential")))
  write_test_osmpbf(path, nodes, ways)

  out_prefix <- tempfile()
  paths <- cpp_extract(path, output = "csv", out_prefix = out_prefix, coords = FALSE, osm = list(progress = FALSE))
  on.exit(unlink(paths), add = TRUE)
  edges_csv <- data.table::fread(paths)
  expect_type(edges_csv$from, "double") # not bit64::integer64
  expect_setequal(c(edges_csv$from, edges_csv$to), c(1, big_id))
})

test_that("cpp_extract on a .pbf resolves oneway correctly and drops non-profile highways", {
  path <- tempfile(fileext = ".pbf")
  on.exit(unlink(path))
  nodes <- data.frame(id = 1:2, lat = c(47.00, 47.01), lon = c(5.00, 5.00))
  ways <- list(list(refs = c(1, 2), tags = list(highway = "footway"))) # not in the default car profile
  write_test_osmpbf(path, nodes, ways)

  res <- cpp_extract(path, output = "data.frame", coords = FALSE, osm = list(progress = FALSE))
  expect_equal(nrow(res$data), 0L)
  expect_null(res$coords)
})

test_that("stream_osmpbf rejects an implausibly large blob header length", {
  path <- make_test_pbf()
  on.exit(unlink(path))
  raw <- readBin(path, "raw", file.info(path)$size)
  bad <- raw
  # First 4 bytes are the BlobHeader length prefix, big-endian -- set it well
  # past the 64 KiB sanity bound (the PBF spec's own cap on BlobHeader size).
  bad[1:4] <- as.raw(c(0x00, 0x10, 0x00, 0x00)) # 1,048,576
  bad_path <- tempfile(fileext = ".pbf")
  on.exit(unlink(bad_path), add = TRUE)
  writeBin(bad, bad_path)
  expect_error(cppRouting:::cpp_osmpbf_counts(bad_path), "implausible blob header length")
})

test_that("stream_osmpbf rejects an implausibly large blob datasize", {
  path <- make_test_pbf()
  on.exit(unlink(path))
  raw <- readBin(path, "raw", file.info(path)$size)
  header_len <- sum(as.integer(raw[1:4]) * c(2^24, 2^16, 2^8, 1))
  header_bytes <- raw[5:(4 + header_len)]
  # BlobHeader body ends with tag(3,0) + varint(datasize) -- its last byte is
  # the (single-byte, since the real blob is tiny) datasize varint. Replace
  # just that value with a properly-encoded varint comfortably past the 64
  # MiB blob-size sanity cap, keeping the preceding tag byte intact, and
  # rewrite the leading big-endian header-length prefix to match the new
  # (longer) header.
  new_datasize_varint <- pbf_varint(100L * 1024L * 1024L)
  new_header <- c(header_bytes[seq_len(length(header_bytes) - 1)], new_datasize_varint)
  new_header_len <- length(new_header)
  len_be <- as.raw(c(bitwAnd(bitwShiftR(new_header_len, 24), 0xFF), bitwAnd(bitwShiftR(new_header_len, 16), 0xFF),
                      bitwAnd(bitwShiftR(new_header_len, 8), 0xFF), bitwAnd(new_header_len, 0xFF)))
  rest <- raw[(5 + header_len):length(raw)]
  bad <- c(len_be, new_header, rest)
  bad_path <- tempfile(fileext = ".pbf")
  on.exit(unlink(bad_path), add = TRUE)
  writeBin(as.raw(bad), bad_path)
  expect_error(cppRouting:::cpp_osmpbf_counts(bad_path), "implausible blob size")
})

test_that("cpp_extract on a .pbf with coords=FALSE drops coords from every output mode", {
  path <- make_test_pbf()
  on.exit(unlink(path))

  res <- cpp_extract(path, output = "data.frame", coords = FALSE, osm = list(progress = FALSE))
  expect_null(res$coords)
  expect_equal(nrow(res$data), 5L)

  g <- cpp_extract(path, coords = FALSE, osm = list(progress = FALSE))
  expect_null(g$coords)

  out_prefix <- tempfile()
  paths <- cpp_extract(path, output = "csv", out_prefix = out_prefix, coords = FALSE, osm = list(progress = FALSE))
  on.exit(unlink(paths), add = TRUE)
  expect_length(paths, 1L)
  expect_true(file.exists(paths))
  expect_false(file.exists(paste0(out_prefix, "_coords.csv")))
})

test_that("cpp_extract on a .pbf stops on geographic coordinates unless coords=FALSE or epsg is set", {
  path <- make_test_pbf()
  on.exit(unlink(path))

  expect_error(cpp_extract(path, output = "data.frame", osm = list(progress = FALSE)), "geographic")

  # coords = FALSE bypasses the issue entirely.
  res <- cpp_extract(path, output = "data.frame", coords = FALSE, osm = list(progress = FALSE))
  expect_null(res$coords)

  # epsg reprojects instead of erroring.
  skip_if_not_installed("sf")
  res2 <- cpp_extract(path, output = "data.frame", epsg = 2154, osm = list(progress = FALSE))
  expect_equal(nrow(res2$coords), 4L)
  expect_true(all(res2$coords$X > 1e5)) # Lambert-93 easting, not longitude degrees

  # csv mode: reprojection happens via a read-back-rewrite of coords.csv.
  out_prefix <- tempfile()
  paths <- cpp_extract(path, output = "csv", out_prefix = out_prefix, epsg = 2154, osm = list(progress = FALSE))
  on.exit(unlink(paths), add = TRUE)
  coords_csv <- data.table::fread(paths[2])
  expect_true(all(coords_csv$X > 1e5))
})

test_that("cpp_extract on a .pbf rejects epsg together with max_memory", {
  out_prefix <- tempfile()
  expect_error(
    cpp_extract("dummy.pbf", output = "csv", out_prefix = out_prefix, epsg = 2154, osm = list(max_memory = 300)),
    "not supported together with max_memory"
  )
  # coords = FALSE sidesteps the conflict (nothing to reproject).
  expect_error(
    cpp_extract("dummy.pbf", output = "csv", out_prefix = out_prefix, epsg = 2154, coords = FALSE, osm = list(max_memory = 300)),
    "File not found" # gets past the epsg/max_memory guard, fails on the (nonexistent) file instead
  )
})

test_that("cpp_extract on a .pbf keep_length keeps haversine length (km) as aux, consistent with dist", {
  path <- make_test_pbf()
  on.exit(unlink(path))

  res <- cpp_extract(path, output = "data.frame", coords = FALSE, keep_length = TRUE, osm = list(progress = FALSE))
  expect_false(is.null(res$aux$length))
  # way A-B-C and way B-D are both "residential" (40 km/h, the default profile speed):
  # dist = 60*length/speed = 1.5*length
  expect_equal(res$data$dist, 1.5 * res$aux$length, tolerance = 1e-6)

  # combined with keep_clazz: both aux columns present
  res2 <- cpp_extract(path, output = "data.frame", coords = FALSE, keep_length = TRUE,
                       osm = list(keep_clazz = TRUE, progress = FALSE))
  expect_setequal(names(res2$aux), c("clazz", "length"))
})

test_that("cpp_extract on a .pbf output='csv' writes the same aux as output='data.frame' to a separate file", {
  path <- make_test_pbf()
  on.exit(unlink(path))

  out_prefix <- tempfile()
  paths <- cpp_extract(path, output = "csv", out_prefix = out_prefix, coords = FALSE, keep_length = TRUE,
                        osm = list(keep_clazz = TRUE, progress = FALSE))
  on.exit(unlink(paths), add = TRUE)
  expect_length(paths, 2L) # edges, aux (no coords)
  aux_csv <- data.table::fread(paths[2])

  res <- cpp_extract(path, output = "data.frame", coords = FALSE, keep_length = TRUE,
                      osm = list(keep_clazz = TRUE, progress = FALSE))
  expect_equal(aux_csv$clazz, res$aux$clazz)
  expect_equal(aux_csv$length, res$aux$length, tolerance = 1e-6)
})
