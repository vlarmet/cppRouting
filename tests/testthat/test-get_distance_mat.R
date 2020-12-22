context("Distance matrices")
library(cppRouting)

edges <- data.frame(from_vertex = c(0, 0, 1, 1, 2, 2, 3, 4, 4),
                    to_vertex = c(1, 3, 2, 4, 4, 5, 1, 3, 5),
                    cost = c(9, 2, 11, 3, 5, 12, 4, 1, 6))
nodes <- unique(c(edges$from_vertex, edges$to_vertex))
directed_graph <- makegraph(edges, directed = TRUE)
non_directed <- makegraph(edges, directed = FALSE)

test_that("directed graph distance matrix", {
    dir_dist <- get_distance_matrix(Graph = directed_graph,
                                    from = nodes, to = nodes, allcores = FALSE)
    output <- matrix(c(0, NA, NA, NA, NA, NA,
                       6, 0, 10, 4, 5, NA,
                       17, 11, 0, 15, 16, NA,
                       2, 4, 6, 0, 1, NA,
                       9, 3, 5, 7, 0, NA,
                       15, 9, 11, 13, 6, 0), 6, 6)
    rownames(output) <- colnames(output) <- c(0:5)
    expect_equal(dir_dist, output)
})

test_that("non directed graph distance matrix", {
    non_dir_dist <- get_distance_matrix(
    Graph = non_directed, from = nodes, to = nodes, allcores = FALSE)
    print(non_dir_dist)
    output <- matrix(c(0, 6, 8, 2, 3, 9,
                       6, 0, 8, 4, 3, 9,
                       8, 8, 0, 6, 5, 11,
                       2, 4, 6, 0, 1, 7,
                       3, 3, 5, 1, 0, 6,
                       9, 9, 11, 7, 6, 0), 6, 6)
    rownames(output) <- colnames(output) <- c(0:5)
    expect_equal(non_dir_dist, output)
})
