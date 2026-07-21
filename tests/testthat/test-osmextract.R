# Confirms the documented recommendation in cpp_extract()'s docs: getting
# the raw cached .pbf path from osmextract (both download_only = TRUE AND
# skip_vectortranslate = TRUE are required -- download_only alone still runs
# the GDAL conversion and returns a .gpkg path instead) and feeding it to
# cpp_extract() like any other .pbf, bypassing osmextract's own sf/gpkg
# conversion entirely.
test_that("the recommended osmextract workflow (raw .pbf path) round-trips", {
  skip_if_not_installed("osmextract")
  skip_if_not_installed("sf")
  test_pbf <- "/home/vincent/cppRouting/pbf/bourgogne-260717.osm.pbf"
  skip_if_not(file.exists(test_pbf), "local test pbf not available")

  path <- osmextract::oe_read(test_pbf, download_only = TRUE, skip_vectortranslate = TRUE, quiet = TRUE)
  expect_true(grepl("\\.osm\\.pbf$", path))

  # epsg required: OSM coordinates are always geographic (WGS84).
  res <- cpp_extract(path, output = "data.frame", epsg = 2154, osm = list(progress = FALSE))
  expect_true(nrow(res$data) > 0)
  expect_true(nrow(res$coords) > 0)
})

test_that("download_only alone (without skip_vectortranslate) returns a .gpkg, not raw .pbf", {
  skip_if_not_installed("osmextract")
  test_pbf <- "/home/vincent/cppRouting/pbf/bourgogne-260717.osm.pbf"
  skip_if_not(file.exists(test_pbf), "local test pbf not available")

  path <- osmextract::oe_read(test_pbf, download_only = TRUE, quiet = TRUE)
  expect_true(grepl("\\.gpkg$", path)) # documents the footgun: this is NOT the raw pbf
})
