# ExternalSorter (src/external_sort.h) is the generic chunked-sort + k-way
# merge building block for the OSM PBF disk pipeline (used when memory can't
# hold the full node/way tables). It knows nothing about OSM -- exercised
# here directly against R's own sort() across chunk sizes small enough to
# force many runs and a real k-way merge, plus edge cases.
test_that("cpp_test_external_sort matches sort() across chunk sizes and edge cases", {
  set.seed(42)
  for (cc in c(1L, 2L, 13L, 100L, 5000L)) {
    # n values at chunk_capacity cc produce ceiling(n/cc) runs, all held
    # open at once during the k-way merge. Windows (MinGW-w64 CRT) caps
    # concurrently open files at ~512 per process, so keep n small enough
    # that even cc = 1 stays well under that; 200 runs is still a genuine
    # k-way merge.
    x <- sample(-1e6:1e6, 200, replace = TRUE) * 1.0
    got <- cppRouting:::cpp_test_external_sort(x, chunk_capacity = cc, temp_dir = tempdir())
    expect_identical(got, sort(x), info = paste("chunk_capacity =", cc))
  }

  expect_identical(cppRouting:::cpp_test_external_sort(numeric(0), 5L, tempdir()), numeric(0))
  expect_identical(cppRouting:::cpp_test_external_sort(42, 5L, tempdir()), 42)
  expect_identical(cppRouting:::cpp_test_external_sort(rep(7, 50), 3L, tempdir()), rep(7, 50))
})

test_that("ExternalSorter cleans up its temp files", {
  before <- list.files(tempdir(), pattern = "test_extsort")
  invisible(cppRouting:::cpp_test_external_sort(sample(1:500, 500), 4L, tempdir()))
  after <- list.files(tempdir(), pattern = "test_extsort")
  expect_equal(before, after)
})
