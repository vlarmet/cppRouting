test_that("basic construction sets nbnode, dict$ref and class", {
  g <- toy_graph()
  expect_s3_class(g, "cppRouting_graph")
  expect_equal(g$nbnode, length(unique(c(toy_edges()$from, toy_edges()$to))))
  expect_setequal(g$dict$ref, as.character(unique(c(toy_edges()$from, toy_edges()$to))))
})

test_that("undirected graphs duplicate edges in both directions", {
  e <- data.frame(from = 1, to = 2, dist = 5)
  g <- makegraph(e, directed = FALSE)
  expect_equal(nrow(g$data), 2)
  expect_equal(get_distance_pair(g, 1, 2), get_distance_pair(g, 2, 1))
})

test_that("df must have exactly 3 columns", {
  expect_error(makegraph(data.frame(a = 1, b = 2)))
  expect_error(makegraph(data.frame(a = 1, b = 2, c = 3, d = 4)))
})

test_that("NAs and negative costs are rejected", {
  expect_error(makegraph(data.frame(from = 1, to = NA, dist = 1)))
  expect_error(makegraph(data.frame(from = 1, to = 2, dist = -1)))
})

test_that("aux as a vector: recycling and length validation", {
  e <- toy_edges()
  g <- makegraph(e, directed = TRUE, aux = 1) # single value recycled
  expect_length(g$attrib$aux, nrow(e))
  expect_true(all(g$attrib$aux == 1))

  expect_error(makegraph(e, directed = TRUE, aux = c(1, 2))) # wrong length
})

test_that("aux as a data.frame: validation rules", {
  e <- toy_edges()

  expect_error(makegraph(e, aux = data.frame())) # 0 columns
  expect_error(makegraph(e, aux = data.frame(x = "a"))) # non-numeric
  expect_error(makegraph(e, aux = setNames(data.frame(1:nrow(e)), ""))) # unnamed
  expect_error(makegraph(e, aux = data.frame(x = 1:2))) # wrong nrow (not 1 or nrow(e))

  # single row recycled to nrow(df)
  g <- makegraph(e, aux = data.frame(x = 1))
  expect_equal(nrow(g$attrib$aux), nrow(e))
})

test_that("negative aux values warn but do not error, for both vector and data.frame aux", {
  e <- toy_edges()
  expect_warning(makegraph(e, aux = c(-1, rep(1, nrow(e) - 1))))
  expect_warning(makegraph(e, aux = data.frame(x = c(-1, rep(1, nrow(e) - 1)))))
})

test_that("capacity/alpha/beta: recycling and validation", {
  e <- toy_edges()
  g <- makegraph(e, capacity = 10, alpha = 0.15, beta = 4)
  expect_length(g$attrib$cap, nrow(e))
  expect_length(g$attrib$alpha, nrow(e))
  expect_length(g$attrib$beta, nrow(e))

  expect_error(makegraph(e, capacity = 0)) # must be strictly positive
  expect_error(makegraph(e, capacity = -1))
  expect_error(makegraph(e, alpha = 0)) # must be != 0
  expect_error(makegraph(e, beta = 0.5)) # must be >= 1
})

test_that("coords: validation and node matching", {
  e <- toy_edges()
  nodes <- unique(c(e$from, e$to))

  good_coords <- data.frame(id = nodes, x = seq_along(nodes), y = seq_along(nodes))
  g <- makegraph(e, coords = good_coords)
  expect_equal(nrow(g$coords), length(nodes))

  expect_error(makegraph(e, coords = data.frame(id = nodes, x = 1))) # wrong ncol
  expect_error(makegraph(e, coords = data.frame(id = c(nodes, nodes[1]), x = 0, y = 0))) # duplicated node
  expect_error(makegraph(e, coords = data.frame(id = nodes[-1], x = 0, y = 0))) # missing node
})
