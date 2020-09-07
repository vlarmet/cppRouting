context("Distance pairs")
library(cppRouting)


test_that("directed graph distance pairs", {
    dir_pairs <- get_distance_pair(Graph = directed_graph,
                                  from = c(0, 0, 0, 0, 0, 1, 2, 3, 4, 5),
                                  to = c(1, 2, 3, 4, 5, 0, 0, 0, 0, 0), allcores = FALSE)
    output <- c(6, 17, 2, 9, 15, NA, NA, NA, NA, NA)
    expect_equal(dir_pairs, output)
})

test_that("non directed graph distance matrix", {
    non_dir_dist <- get_distance_pair(Graph = non_directed,
                                      from = c(0, 0, 0, 0, 0, 1, 2, 3, 4, 5),
                                      to = c(1, 2, 3, 4, 5, 0, 0, 0, 0, 0), allcores = FALSE)
    print(non_dir_dist)
    output <- c(6, 8, 2, 3, 9, 6, 8, 2, 3, 9)
    expect_equal(non_dir_dist, output)
})
