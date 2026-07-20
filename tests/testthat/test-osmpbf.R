test_that("cpp_osmpbf_counts reads back the exact synthetic fixture", {
  path <- make_test_pbf()
  on.exit(unlink(path))
  counts <- cppRouting:::cpp_osmpbf_counts(path)
  expect_equal(counts$n_nodes, 4)
  expect_equal(counts$n_ways, 2)
  expect_equal(counts$n_tagged_ways, 2)
})

test_that("cpp_extract on a .pbf splits a way at a mid-way junction (node 2)", {
  path <- make_test_pbf()
  on.exit(unlink(path))

  res <- cpp_extract(path, output = "data.frame", progress = FALSE)

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

test_that("cpp_extract on a .pbf resolves oneway correctly and drops non-profile highways", {
  path <- tempfile(fileext = ".pbf")
  on.exit(unlink(path))
  nodes <- data.frame(id = 1:2, lat = c(47.00, 47.01), lon = c(5.00, 5.00))
  ways <- list(list(refs = c(1, 2), tags = list(highway = "footway"))) # not in the default car profile
  write_test_osmpbf(path, nodes, ways)

  res <- cpp_extract(path, output = "data.frame", progress = FALSE)
  expect_equal(nrow(res$data), 0L)
  expect_equal(nrow(res$coords), 0L)
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
