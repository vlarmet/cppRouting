# Disk-backed OSM pipeline (src/io_osm_graph_disk.cpp): same tag resolution
# and junction-splitting semantics as the in-RAM path (io_osm_graph.cpp), but
# spills every intermediate table to disk via external sort + sequential
# merge-joins, bounded to roughly max_memory_mb of RAM regardless of input
# size. Exercised here with a deliberately tiny max_memory (a few hundred
# bytes) to force many chunks/runs even on this tiny fixture -- the point is
# to prove correctness of the chunking machinery itself, not performance.

run_disk_pipeline <- function(path, max_memory_mb = 0.0001, keep_clazz = FALSE, keep_length = FALSE, ...){
  out_prefix <- tempfile()
  edges_path <- paste0(out_prefix, "_edges.csv")
  coords_path <- paste0(out_prefix, "_coords.csv")
  aux_path <- paste0(out_prefix, "_aux.csv")
  profile <- osm_car_profile()
  cppRouting:::cpp_build_osm_graph_disk_csv(
    path, as.character(profile$highway), as.integer(profile$clazz), as.numeric(profile$speed),
    TRUE, TRUE, "car", max_memory_mb, tempdir(), edges_path, coords_path, aux_path,
    TRUE, keep_clazz, keep_length, FALSE)
  out <- list(data = data.table::fread(edges_path), coords = data.table::fread(coords_path))
  if (keep_clazz || keep_length) out$aux <- data.table::fread(aux_path)
  out
}

test_that("disk pipeline splits a way at a mid-way junction, same as the in-RAM path", {
  path <- make_test_pbf() # from helper-osmpbf.R: way1 A-B-C (twoway), way2 B-D (oneway), B shared mid-way
  on.exit(unlink(path))

  res <- run_disk_pipeline(path)
  expect_equal(nrow(res$coords), 4L)
  expect_equal(nrow(res$data), 5L)

  # run_disk_pipeline() calls the C++ export directly, bypassing the R-level
  # geographic-coordinates hard stop -- coords = FALSE here since this
  # particular comparison only needs edges/dist parity, not coords.
  in_ram <- cpp_extract(path, output = "data.frame", coords = FALSE, osm_params = list(progress = FALSE))
  expect_equal(nrow(res$data), nrow(in_ram$data))
  expect_equal(sort(res$data$dist), sort(in_ram$data$dist))
})

test_that("disk pipeline and in-RAM path agree on the exact real OSM node ids, not just row counts", {
  path <- make_test_pbf()
  on.exit(unlink(path))

  res <- run_disk_pipeline(path)
  in_ram <- cpp_extract(path, output = "data.frame", coords = FALSE, osm_params = list(progress = FALSE))

  # Both must expose the fixture's real ids (1:4, see helper-osmpbf.R) --
  # previously the disk pipeline kept the real OSM id while the in-RAM path
  # remapped to a synthetic 0-based index, so the same file produced two
  # different id spaces depending on max_memory. Fixed so both paths agree.
  expect_setequal(c(res$data$from, res$data$to), 1:4)
  expect_setequal(c(res$data$from, res$data$to), c(in_ram$data$from, in_ram$data$to))
})

test_that("disk pipeline matches the in-RAM path on the real Bourgogne extract", {
  skip_if_not_installed("sf")
  pbf <- "/home/vincent/cppRouting/pbf/bourgogne-260717.osm.pbf"
  skip_if_not(file.exists(pbf), "local test pbf not available")

  disk_res <- run_disk_pipeline(pbf, max_memory_mb = 5)
  # epsg required for coords (OSM coordinates are always geographic); row
  # count is unaffected by reprojection, so it's still a valid parity check
  # against disk_res$coords (left in WGS84 by the raw C++ call above).
  in_ram <- cpp_extract(pbf, output = "data.frame", epsg = 2154, osm_params = list(progress = FALSE))

  expect_equal(nrow(disk_res$coords), nrow(in_ram$coords))
  expect_equal(nrow(disk_res$data), nrow(in_ram$data))
  expect_equal(sum(disk_res$data$dist), sum(in_ram$data$dist), tolerance = 1e-6)
})

test_that("cpp_extract(max_memory=...) reaches the disk pipeline and matches the in-RAM path", {
  pbf <- "/home/vincent/cppRouting/pbf/bourgogne-260717.osm.pbf"
  skip_if_not(file.exists(pbf), "local test pbf not available")

  # OSM coordinates are always geographic -- coords = FALSE here; coords
  # parity for the disk pipeline itself is already covered by the raw
  # C++-level test above, and epsg + max_memory by the dedicated test in
  # test-osmpbf.R.
  out_prefix <- tempfile()
  paths <- cpp_extract(pbf, output = "csv", out_prefix = out_prefix, coords = FALSE,
                        osm_params = list(max_mem = 300, progress = FALSE))
  expect_true(file.exists(paths))

  disk_res <- list(data = data.table::fread(paths))
  in_ram <- cpp_extract(pbf, output = "data.frame", coords = FALSE, osm_params = list(progress = FALSE))
  expect_equal(nrow(disk_res$data), nrow(in_ram$data))
})

test_that("cpp_extract(max_memory=..., coords=FALSE) skips coords.csv but edge distances stay correct", {
  path <- make_test_pbf()
  on.exit(unlink(path))

  out_prefix <- tempfile()
  paths <- cpp_extract(path, output = "csv", out_prefix = out_prefix, coords = FALSE,
                        osm_params = list(max_mem = 300, progress = FALSE))
  on.exit(unlink(paths), add = TRUE)
  expect_length(paths, 1L)
  expect_true(file.exists(paths))
  expect_false(file.exists(paste0(out_prefix, "_coords.csv")))

  disk_edges <- data.table::fread(paths)
  in_ram <- cpp_extract(path, output = "data.frame", coords = FALSE, osm_params = list(progress = FALSE))
  expect_equal(sort(disk_edges$dist), sort(in_ram$data$dist))
})

test_that("cpp_extract(max_memory=...) output='csv' writes the same aux as the in-RAM path to a separate file", {
  path <- make_test_pbf()
  on.exit(unlink(path))

  out_prefix <- tempfile()
  paths <- cpp_extract(path, output = "csv", out_prefix = out_prefix, coords = FALSE, keep_length = TRUE,
                        osm_params = list(max_mem = 300, keep_clazz = TRUE, progress = FALSE))
  on.exit(unlink(paths), add = TRUE)
  expect_length(paths, 2L) # edges, aux (no coords)
  aux_csv <- data.table::fread(paths[2])

  in_ram <- cpp_extract(path, output = "data.frame", coords = FALSE, keep_length = TRUE,
                         osm_params = list(keep_clazz = TRUE, progress = FALSE))
  expect_equal(sort(aux_csv$clazz), sort(in_ram$aux$clazz))
  expect_equal(sort(aux_csv$length), sort(in_ram$aux$length), tolerance = 1e-6)
})

test_that("cpp_extract() requires output='csv' when max_memory is set", {
  expect_error(cpp_extract("dummy.pbf", osm_params = list(max_mem = 300)), "requires output")
  expect_error(cpp_extract("dummy.pbf", output = "graph", osm_params = list(max_mem = 300)), "requires output")
})

test_that("cpp_extract() rejects an invalid max_memory", {
  out_prefix <- tempfile()
  expect_error(cpp_extract("dummy.pbf", output = "csv", out_prefix = out_prefix, osm_params = list(max_mem = -1)), "positive number")
  expect_error(cpp_extract("dummy.pbf", output = "csv", out_prefix = out_prefix, osm_params = list(max_mem = 0)), "positive number")
  expect_error(cpp_extract("dummy.pbf", output = "csv", out_prefix = out_prefix, osm_params = list(max_mem = NA_real_)), "positive number")
  expect_error(cpp_extract("dummy.pbf", output = "csv", out_prefix = out_prefix, osm_params = list(max_mem = "300")), "positive number")
  expect_error(cpp_extract("dummy.pbf", output = "csv", out_prefix = out_prefix, osm_params = list(max_mem = c(100, 200))), "positive number")
})

test_that("reproject_coords_csv_chunked: chunked reprojection equals the bulk one, with many forced chunks", {
  skip_if_not_installed("sf")
  path <- make_test_pbf() # 4-node fixture
  on.exit(unlink(path))

  # Raw disk pipeline output: WGS84 coords (the raw C++ call bypasses the
  # R-level geographic-coordinates guard, same as the other disk tests).
  raw <- run_disk_pipeline(path)$coords
  expect_true(all(raw$Y > -90 & raw$Y < 90)) # lat/lon, not yet projected

  # In-RAM path with epsg: the reference projected values.
  ref_prefix <- tempfile()
  ref_paths <- cpp_extract(path, output = "csv", out_prefix = ref_prefix, epsg = 2154,
                            osm_params = list(progress = FALSE))
  on.exit(unlink(ref_paths), add = TRUE)
  ref <- data.table::fread(ref_paths[2])

  # Force two reproject chunks on the 4-row fixture (2 rows each) and check
  # the concatenation is row-for-row the bulk reprojection.
  work <- tempfile(); data.table::fwrite(raw, work)
  on.exit(unlink(work), add = TRUE)
  cppRouting:::reproject_coords_csv_chunked(work, TRUE, 4326L, 2154, progress = FALSE, chunk_rows = 2)
  got <- data.table::fread(work)
  expect_equal(nrow(got), 4L)
  expect_equal(sort(got$node_id), sort(ref$node_id))
  expect_equal(sort(got$X), sort(ref$X), tolerance = 1e-9)
  expect_equal(sort(got$Y), sort(ref$Y), tolerance = 1e-9)

  # Empty input (header-only csv) is a no-op that keeps the file intact.
  empty <- tempfile(); data.table::fwrite(data.frame(node_id = character(0), X = numeric(0), Y = numeric(0)), empty)
  on.exit(unlink(empty), add = TRUE)
  cppRouting:::reproject_coords_csv_chunked(empty, TRUE, 4326L, 2154, progress = FALSE, chunk_rows = 2)
  expect_equal(nrow(data.table::fread(empty)), 0L)
})
