# Disk-backed OSM pipeline (src/io_osm_graph_disk.cpp): same tag resolution
# and junction-splitting semantics as the in-RAM path (io_osm_graph.cpp), but
# spills every intermediate table to disk via external sort + sequential
# merge-joins, bounded to roughly max_memory_mb of RAM regardless of input
# size. Exercised here with a deliberately tiny max_memory (a few hundred
# bytes) to force many chunks/runs even on this tiny fixture -- the point is
# to prove correctness of the chunking machinery itself, not performance.

run_disk_pipeline <- function(path, max_memory_mb = 0.0001, ...){
  out_prefix <- tempfile()
  edges_path <- paste0(out_prefix, "_edges.csv")
  coords_path <- paste0(out_prefix, "_coords.csv")
  profile <- osm_car_profile()
  cppRouting:::cpp_build_osm_graph_disk_csv(
    path, as.character(profile$highway), as.integer(profile$clazz), as.numeric(profile$speed),
    TRUE, TRUE, max_memory_mb, tempdir(), edges_path, coords_path, TRUE, FALSE)
  list(data = data.table::fread(edges_path), coords = data.table::fread(coords_path))
}

test_that("disk pipeline splits a way at a mid-way junction, same as the in-RAM path", {
  path <- make_test_pbf() # from helper-osmpbf.R: way1 A-B-C (twoway), way2 B-D (oneway), B shared mid-way
  on.exit(unlink(path))

  res <- run_disk_pipeline(path)
  expect_equal(nrow(res$coords), 4L)
  expect_equal(nrow(res$data), 5L)

  in_ram <- cpp_extract(path, output = "data.frame", progress = FALSE)
  expect_equal(nrow(res$data), nrow(in_ram$data))
  expect_equal(sort(res$data$dist), sort(in_ram$data$dist))
})

test_that("disk pipeline matches the in-RAM path on the real Bourgogne extract", {
  pbf <- "/home/vincent/cppRouting/pbf/bourgogne-260717.osm.pbf"
  skip_if_not(file.exists(pbf), "local test pbf not available")

  disk_res <- run_disk_pipeline(pbf, max_memory_mb = 5)
  in_ram <- cpp_extract(pbf, output = "data.frame", progress = FALSE)

  expect_equal(nrow(disk_res$coords), nrow(in_ram$coords))
  expect_equal(nrow(disk_res$data), nrow(in_ram$data))
  expect_equal(sum(disk_res$data$dist), sum(in_ram$data$dist), tolerance = 1e-6)
})

test_that("cpp_extract(max_memory=...) reaches the disk pipeline and matches the in-RAM path", {
  pbf <- "/home/vincent/cppRouting/pbf/bourgogne-260717.osm.pbf"
  skip_if_not(file.exists(pbf), "local test pbf not available")

  out_prefix <- tempfile()
  paths <- cpp_extract(pbf, max_memory = 300, output = "csv", out_prefix = out_prefix, progress = FALSE)
  expect_true(all(file.exists(paths)))

  disk_res <- list(data = data.table::fread(paths[1]), coords = data.table::fread(paths[2]))
  in_ram <- cpp_extract(pbf, output = "data.frame", progress = FALSE)
  expect_equal(nrow(disk_res$coords), nrow(in_ram$coords))
  expect_equal(nrow(disk_res$data), nrow(in_ram$data))
})

test_that("cpp_extract(max_memory=..., coords=FALSE) skips coords.csv but edge distances stay correct", {
  path <- make_test_pbf()
  on.exit(unlink(path))

  out_prefix <- tempfile()
  paths <- cpp_extract(path, max_memory = 300, output = "csv", out_prefix = out_prefix, coords = FALSE, progress = FALSE)
  on.exit(unlink(paths), add = TRUE)
  expect_length(paths, 1L)
  expect_true(file.exists(paths))
  expect_false(file.exists(paste0(out_prefix, "_coords.csv")))

  disk_edges <- data.table::fread(paths)
  in_ram <- cpp_extract(path, output = "data.frame", progress = FALSE)
  expect_equal(sort(disk_edges$dist), sort(in_ram$data$dist))
})

test_that("cpp_extract() requires output='csv' when max_memory is set", {
  expect_error(cpp_extract("dummy.pbf", max_memory = 300), "requires output")
  expect_error(cpp_extract("dummy.pbf", max_memory = 300, output = "graph"), "requires output")
})

test_that("cpp_extract() rejects an invalid max_memory", {
  out_prefix <- tempfile()
  expect_error(cpp_extract("dummy.pbf", max_memory = -1, output = "csv", out_prefix = out_prefix), "positive number")
  expect_error(cpp_extract("dummy.pbf", max_memory = 0, output = "csv", out_prefix = out_prefix), "positive number")
  expect_error(cpp_extract("dummy.pbf", max_memory = NA_real_, output = "csv", out_prefix = out_prefix), "positive number")
  expect_error(cpp_extract("dummy.pbf", max_memory = "300", output = "csv", out_prefix = out_prefix), "positive number")
  expect_error(cpp_extract("dummy.pbf", max_memory = c(100, 200), output = "csv", out_prefix = out_prefix), "positive number")
})
