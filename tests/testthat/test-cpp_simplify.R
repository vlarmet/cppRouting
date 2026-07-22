# Uses the exact toy graph from README.md's "Network simplification" section,
# where the expected behavior at each step is documented by the package
# author - verified empirically to match before writing these assertions.
simplify_demo_edges <- function() {
  data.frame(
    from = c("a", "b", "c", "d", "d", "e", "e", "e", "f", "f", "g", "h", "h", "h",
             "i", "j", "k", "k", "k", "l", "l", "l", "m", "m", "m", "n", "n", "o", "p", "q", "r"),
    to   = c("b", "c", "d", "e", "k", "f", "d", "h", "g", "e", "f", "e", "i", "k",
             "j", "i", "h", "d", "l", "k", "m", "n", "n", "o", "l", "l", "m", "m", "r", "p", "q"),
    dist = rep(1, 31)
  )
}

test_that("one pass contracts non-intersection nodes (b, c, f, n)", {
  g <- makegraph(simplify_demo_edges(), directed = TRUE)
  s <- cpp_simplify(g, rm_loop = FALSE)
  nodes <- unique(c(to_df(s)$from, to_df(s)$to))
  expect_false(any(c("b", "c", "f", "n") %in% nodes))
  expect_true(all(c("a", "d", "e", "g", "h", "i", "j", "k", "l", "m", "o", "p", "q", "r") %in% nodes))
})

test_that("iterate = TRUE removes nodes that only become non-intersection after the first pass (l, m)", {
  g <- makegraph(simplify_demo_edges(), directed = TRUE)
  s_once <- cpp_simplify(g, rm_loop = FALSE, iterate = FALSE)
  s_iter <- cpp_simplify(g, rm_loop = FALSE, iterate = TRUE)

  expect_true(all(c("l", "m") %in% unique(c(to_df(s_once)$from, to_df(s_once)$to))))
  expect_false(any(c("l", "m") %in% unique(c(to_df(s_iter)$from, to_df(s_iter)$to))))
  # the collapsed k->l->m->o chain (cost 1+1+1) becomes a single k->o edge
  merged <- to_df(s_iter)
  expect_equal(merged$dist[merged$from == "k" & merged$to == "o"], 3)
})

test_that("rm_loop removes a fully isolated cycle instead of leaving a self-loop", {
  # p->r->q->p is an isolated 3-cycle disconnected from the rest of the graph
  g <- makegraph(simplify_demo_edges(), directed = TRUE)
  s_keep_loop <- cpp_simplify(g, rm_loop = FALSE, iterate = TRUE)
  s_rm_loop <- cpp_simplify(g, rm_loop = TRUE, iterate = TRUE)

  expect_true(all(c("p", "q", "r") %in% unique(c(to_df(s_keep_loop)$from, to_df(s_keep_loop)$to))))
  expect_false(any(c("p", "q", "r") %in% unique(c(to_df(s_rm_loop)$from, to_df(s_rm_loop)$to))))
})

test_that("keep protects specified nodes from being contracted", {
  g <- makegraph(simplify_demo_edges(), directed = TRUE)
  s <- cpp_simplify(g, rm_loop = FALSE, keep = "b")
  nodes <- unique(c(to_df(s)$from, to_df(s)$to))
  expect_true("b" %in% nodes)
})

test_that("cpp_simplify() result is classed and usable by downstream functions", {
  g <- grid_graph(5, 5, directed = TRUE)
  s <- cpp_simplify(g, iterate = TRUE)
  expect_s3_class(s, "cppRouting_simplified")
  expect_s3_class(s, "cppRouting_graph")
  expect_no_error(get_distance_pair(s, "1", "25"))
  expect_no_error(cpp_contract(s, silent = TRUE))
})
