# Shared fixtures for the test suite. Kept small and deterministic (no real
# road data - that lives in data_readme/ and is used by the bench_*.R
# performance scripts instead, not by unit tests).

# A hand-verifiable directed graph, small enough that shortest paths/costs
# can be checked by inspection. Includes:
# - a parallel edge pair (1->2) with different costs, to exercise
#   "cheapest parallel edge wins" (a real bug class hit repeatedly while
#   building aggregate_aux)
# - an exact tie: two paths of equal total weight but different edge counts
#   (1->4 via 1-2-4 and 1-3-4, both cost 5, hop count 2 either way here -
#   kept simple; ties-with-different-hop-count are exercised in the bigger
#   fixture below since they need more structure)
toy_edges <- function() {
  data.frame(
    from = c(1, 1, 1, 2, 2, 3, 4, 4, 5),
    to   = c(2, 2, 3, 4, 5, 4, 5, 6, 6),
    dist = c(2, 5, 3, 3, 8, 2, 1, 6, 1),
    stringsAsFactors = FALSE
  )
}

toy_aux <- function() {
  e <- toy_edges()
  data.frame(hop = rep(1, nrow(e)), dist2 = e$dist * 2)
}

toy_graph <- function(directed = TRUE, aux = NULL) {
  cppRouting::makegraph(toy_edges(), directed = directed, aux = aux)
}

# A small grid graph: realistic/near-planar topology (unlike a fully random
# graph, which is an adversarial worst case for both CH and CCH elimination
# ordering - see min_degree_order fill-in blowup found during development).
# Big enough to exercise contraction/CCH/simplification meaningfully, small
# enough to run fast. Nodes are numbered row-major, 1-based.
grid_edges <- function(nrow = 6, ncol = 6, directed = TRUE, seed = 1) {
  set.seed(seed)
  node_id <- function(r, c) (r - 1) * ncol + c
  from <- integer(0); to <- integer(0)
  for (r in seq_len(nrow)) {
    for (c in seq_len(ncol)) {
      if (c < ncol) { from <- c(from, node_id(r, c)); to <- c(to, node_id(r, c + 1)) }
      if (r < nrow) { from <- c(from, node_id(r, c)); to <- c(to, node_id(r + 1, c)) }
    }
  }
  weight <- round(runif(length(from), 1, 10), 3)
  data.frame(from = from, to = to, weight = weight)
}

grid_coords <- function(nrow = 6, ncol = 6) {
  node_id <- function(r, c) (r - 1) * ncol + c
  ids <- integer(0); x <- numeric(0); y <- numeric(0)
  for (r in seq_len(nrow)) for (c in seq_len(ncol)) {
    ids <- c(ids, node_id(r, c)); x <- c(x, c); y <- c(y, r)
  }
  data.frame(id = ids, x = x, y = y)
}

# Grid graph with two aux columns (a constant "hop" column, and one that
# mirrors the main weight), coordinates, and capacity/alpha/beta so it can
# be reused for aggregate_aux, A*/NBA and traffic-assignment tests alike.
grid_graph <- function(nrow = 6, ncol = 6, directed = TRUE, seed = 1,
                       with_coords = FALSE, with_aux = FALSE,
                       with_traffic_attrs = FALSE) {
  e <- grid_edges(nrow, ncol, directed, seed)
  args <- list(df = e, directed = directed)
  if (with_coords) args$coords <- grid_coords(nrow, ncol)
  if (with_aux) args$aux <- data.frame(hop = rep(1, nrow(e)), dist2 = e$weight)
  if (with_traffic_attrs) {
    args$capacity <- round(runif(nrow(e), 50, 200), 1)
    args$alpha <- rep(0.15, nrow(e))
    args$beta <- rep(4, nrow(e))
  }
  do.call(cppRouting::makegraph, args)
}

# Duplicated edge with a strictly higher cost than its twin, to make sure
# tests exercise "cheapest parallel edge wins" rather than an arbitrary one.
with_parallel_edge <- function(edges, from, to, extra_cost = 5) {
  base <- edges[edges$from == from & edges$to == to, ]
  stopifnot(nrow(base) >= 1)
  dup <- base[1, ]
  dup[[3]] <- dup[[3]] + extra_cost
  rbind(edges, dup)
}
