context("Distance matrices")
library(cppRouting)

test_that("make directed graph", {
  x <- c(0, 0, 1, 1, 2, 2, 3, 4, 4)
  y <- c(1, 3, 2, 4, 4, 5, 1, 3, 5)
  cost <- c(9, 2, 11, 3, 5, 12, 4, 1, 6)
  edges <- data.frame(from_vertex = LETTERS[x + 1],
                      to_vertex = LETTERS[y + 1],
                      cost = cost)
  directed_graph <- makegraph(edges, directed = TRUE)
  expect_equal(directed_graph$data$from, x)
  expect_equal(directed_graph$data$to, y)
  expect_equal(directed_graph$data$dist, cost)
})

test_that("make undirected graph", {
  x <- c(0, 0, 1, 1, 2, 2, 3, 4, 4)
  y <- c(1, 3, 2, 4, 4, 5, 1, 3, 5)
  cost <- c(9, 2, 11, 3, 5, 12, 4, 1, 6)
  edges <- data.frame(from_vertex = LETTERS[x + 1],
                      to_vertex = LETTERS[y + 1],
                      cost = cost)
  udg <- makegraph(edges, directed = FALSE)
  expect_equal(udg$data[udg$data$from == 1 & udg$data$to == 0,]$dist, 9)
})

test_that("make directed graph with coordinates", {
  x <- c(0, 0, 1, 1, 2, 2, 3, 4, 4)
  y <- c(1, 3, 2, 4, 4, 5, 1, 3, 5)
  cost <- c(9, 2, 11, 3, 5, 12, 4, 1, 6)
  edges <- data.frame(from_vertex = LETTERS[x + 1],
                    to_vertex = LETTERS[y + 1],
                    cost = cost)
  coord <- data.frame(node = LETTERS[c(1:6)],
                      X = c(2, 2, 2, 0, 0, 0),
                      Y = c(0, 2, 2, 0, 2, 4))
  directed_graph <- makegraph(edges, directed = TRUE, coords = coord)
  print(directed_graph$nodes)
  expect_equal(directed_graph$coords, coord)
  expect_equal(directed_graph$nbnode, 6)
  expect_equal(directed_graph$dict, data.frame(ref = LETTERS[c(1:6)], id = c(0:5)))

})



test_that("coords not unique", {
  x <- c(0, 0, 1, 1, 2, 2, 3, 4, 4)
  y <- c(1, 3, 2, 4, 4, 5, 1, 3, 5)
  cost <- c(9, 2, 11, 3, 5, 12, 4, 1, 6)
  edges <- data.frame(from_vertex = LETTERS[x + 1],
                    to_vertex = LETTERS[y + 1],
                    cost = cost
                    )
  coord <- data.frame(node = LETTERS[c(1:5, 5)],
                      X = c(2, 2, 2, 0, 0, 0),
                      Y = c(0, 2, 2, 0, 2, 4))

  expect_error(makegraph(edges, directed = TRUE, coords = coord),
               "Nodes should be unique in the coordinates data frame")
  })


test_that("missing nodes", {
  x <- c(0, 0, 1, 1, 2, 2, 3, 4, 4)
  y <- c(1, 3, 2, 4, 4, 5, 1, 3, 5)
  cost <- c(9, 2, 11, 3, 5, 12, 4, 1, 6)
  edges <- data.frame(from_vertex = LETTERS[x + 1],
                    to_vertex = LETTERS[y + 1],
                    cost = cost
                    )
  coord <- data.frame(node = LETTERS[c(1:5)],
                      X = c(2, 2, 2, 0, 0),
                      Y = c(0, 2, 2, 0, 2))

  expect_error(makegraph(edges, directed = TRUE, coords = coord),
               "Some nodes are missing in coordinates data")
})


test_that("NAs in coords", {
  x <- c(0, 0, 1, 1, 2, 2, 3, 4, 4)
  y <- c(1, 3, 2, 4, 4, 5, 1, 3, 5)
  cost <- c(9, 2, 11, 3, 5, 12, 4, 1, 6)
  edges <- data.frame(from_vertex = LETTERS[x + 1],
                    to_vertex = LETTERS[y + 1],
                    cost = cost
                    )
  coord <- data.frame(node = LETTERS[c(1:6)],
                      X = c(2, 2, 2, 0, 0, 0),
                      Y = c(0, NA, 2, 0, 2, 4))

  expect_error(makegraph(edges, directed = TRUE, coords = coord),
               "NAs are not allowed in coordinates")
  coord <- data.frame(node = LETTERS[c(1:6)],
                      X = c(NA, 2, 2, 0, 0, 0),
                      Y = c(0, 2, 2, 0, 2, 4))
  expect_error(makegraph(edges, directed = TRUE, coords = coord),
               "NAs are not allowed in coordinates")
})

test_that("Wrong size data object supplied", {
  edges <- data.frame(w = c(1), x = c(2), y = c(3), z = c(4))
  expect_error(makegraph(edges), "Data should have 3 columns")
})

test_that("NAs rejected in graph", {
  edges <- data.frame(x = c(NA), y = c(2), z = c(3))
  expect_error(makegraph(edges), "NAs are not allowed in the graph")
  edges <- data.frame(x = c(1), y = c(NA), z = c(3))
  expect_error(makegraph(edges), "NAs are not allowed in the graph")
})

test_that("Negative costs rejected", {
  edges <- data.frame(x = c(1), y = c(2), z = c(-3))
  expect_error(makegraph(edges), "Negative cost is not allowed")
})
