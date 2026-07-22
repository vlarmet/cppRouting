test_that("to_df() returns only from/to/dist by default, aux/capacity/alpha/beta all FALSE", {
  g <- toy_graph()
  d <- to_df(g)
  expect_named(d, c("from", "to", "dist"))
})

test_that("to_df(aux = TRUE) appends a plain vector aux as a column named 'aux'", {
  e <- toy_edges()
  g <- makegraph(e, directed = TRUE, aux = e$dist * 10)
  d <- to_df(g, aux = TRUE)
  expect_true("aux" %in% names(d))
  expect_equal(d$aux, g$attrib$aux)

  # not requested: no aux column even though the graph has one
  expect_false("aux" %in% names(to_df(g)))
})

test_that("to_df(aux = TRUE) appends a data.frame aux keeping its own column names", {
  e <- toy_edges()
  g <- makegraph(e, directed = TRUE, aux = toy_aux())
  d <- to_df(g, aux = TRUE)
  expect_true(all(c("hop", "dist2") %in% names(d)))
  expect_equal(d$hop, g$attrib$aux$hop)
  expect_equal(d$dist2, g$attrib$aux$dist2)
})

test_that("to_df(capacity/alpha/beta = TRUE) append their respective columns", {
  e <- toy_edges()
  g <- makegraph(e, directed = TRUE, capacity = 10, alpha = 0.15, beta = 4)
  d <- to_df(g, capacity = TRUE, alpha = TRUE, beta = TRUE)
  expect_equal(d$capacity, g$attrib$cap)
  expect_equal(d$alpha, g$attrib$alpha)
  expect_equal(d$beta, g$attrib$beta)
})

test_that("to_df() silently ignores aux/capacity/alpha/beta requested but absent from the graph", {
  g <- toy_graph() # no aux/capacity/alpha/beta at all
  d <- to_df(g, aux = TRUE, capacity = TRUE, alpha = TRUE, beta = TRUE)
  expect_named(d, c("from", "to", "dist")) # no extra columns, no error
})
