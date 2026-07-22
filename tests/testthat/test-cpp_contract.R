test_that("customizable = FALSE (default) produces a classical CH graph", {
  g <- toy_graph()
  ch <- cpp_contract(g, silent = TRUE)
  expect_s3_class(ch, "cppRouting_contracted")
  expect_false(inherits(ch, "cppRouting_cch"))
})

test_that("customizable = TRUE produces a CCH topology object", {
  g <- toy_graph()
  cch <- cpp_contract(g, customizable = TRUE, silent = TRUE)
  expect_s3_class(cch, "cppRouting_cch")
  expect_false(inherits(cch, "cppRouting_contracted"))
})

test_that("only makegraph()/cpp_simplify() output is accepted", {
  expect_error(cpp_contract(list(a = 1)))
  expect_error(cpp_contract(data.frame(from = 1, to = 2, dist = 3)))
})

test_that("a custom order (customizable = TRUE) still produces correct distances", {
  g <- grid_graph(5, 5, directed = TRUE)
  set.seed(31)
  custom_order <- sample(g$dict$ref)
  cch_default <- cpp_customize(cpp_contract(g, customizable = TRUE, silent = TRUE))
  cch_ordered <- cpp_customize(cpp_contract(g, customizable = TRUE, order = custom_order, silent = TRUE))

  from <- as.character(1:5)
  to <- as.character(21:25)
  expect_equal(get_distance_pair(cch_default, from, to), get_distance_pair(cch_ordered, from, to), tolerance = 1e-8)
})

test_that("order is validated: wrong length, duplicates and unknown nodes are rejected", {
  g <- grid_graph(4, 4, directed = TRUE)
  nodes <- g$dict$ref

  expect_error(cpp_contract(g, customizable = TRUE, order = nodes[1:3])) # wrong length

  dup_order <- nodes
  dup_order[length(dup_order)] <- dup_order[1] # nodes[1] twice, last node missing
  expect_error(cpp_contract(g, customizable = TRUE, order = dup_order))

  unknown_order <- nodes
  unknown_order[1] <- "not_a_node"
  expect_error(cpp_contract(g, customizable = TRUE, order = unknown_order))
})

test_that("silent = TRUE suppresses progress output for both contraction flavors", {
  g <- grid_graph(5, 5, directed = TRUE)
  expect_silent(cpp_contract(g, silent = TRUE))
  expect_silent(cpp_contract(g, customizable = TRUE, silent = TRUE))
})
