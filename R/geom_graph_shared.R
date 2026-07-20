# Internal helpers shared by graph_from_shapefile(), graph_from_geopackage()
# and graph_from_sf() (Module A: "geometry endpoint graph builder", see the
# graph-import plan). Each format-specific reader only has to produce a flat
# per-segment table (x1,y1,x2,y2,length,attrs) -- everything from cost/
# direction resolution down to the final graph/data.frame/csv output is
# common and lives here.

resolve_geom_cost <- function(attrs, geom_length, cost_field, speed_field, length_field){
  if (!is.null(cost_field)) {
    if (!cost_field %in% names(attrs)) stop("cost_field '", cost_field, "' not found in attributes")
    cost <- as.numeric(attrs[[cost_field]])
  } else if (!is.null(speed_field)) {
    if (!speed_field %in% names(attrs)) stop("speed_field '", speed_field, "' not found in attributes")
    length_val <- if (!is.null(length_field)) {
      if (!length_field %in% names(attrs)) stop("length_field '", length_field, "' not found in attributes")
      as.numeric(attrs[[length_field]])
    } else geom_length
    cost <- 60 * length_val / as.numeric(attrs[[speed_field]])
  } else {
    cost <- geom_length
  }
  if (any(is.na(cost)) || any(cost < 0)) stop("Resolved cost contains NA or negative values")
  cost
}

resolve_geom_direction <- function(attrs, direction_field, direction_values, directed, n){
  direction_code <- rep(if (directed) 1L else 0L, n)
  if (is.null(direction_field)) return(direction_code)

  if (!direction_field %in% names(attrs)) stop("direction_field '", direction_field, "' not found in attributes")
  if (is.null(direction_values)) stop("direction_values must be provided together with direction_field")

  dvals <- as.character(attrs[[direction_field]])
  direction_code <- rep(0L, n)
  if (!is.null(direction_values$forward)) direction_code[dvals %in% direction_values$forward] <- 1L
  if (!is.null(direction_values$backward)) direction_code[dvals %in% direction_values$backward] <- 2L

  known <- unlist(direction_values[c("forward", "backward", "twoway")], use.names = FALSE)
  unmatched <- !dvals %in% known
  if (any(unmatched)) {
    warning(sum(unmatched), " edge(s) have a direction_field value not listed in direction_values; ",
            "treated as twoway (", paste(unique(dvals[unmatched]), collapse = ", "), ")")
  }
  direction_code
}

warn_if_geographic_without_length <- function(is_geographic, length_field, cost_field){
  if (isTRUE(is_geographic) && is.null(length_field) && is.null(cost_field)) {
    warning("Coordinates appear to be geographic (degrees), not projected. ",
            "Geometric length computed from these coordinates will not be a reliable planar distance; ",
            "consider reprojecting first or supplying length_field/cost_field directly.")
  }
}

# Dispatches to either the CSV streaming writer (edges/coords written
# straight to disk as they're produced, aux never carried -- same
# restriction as the OSM CSV path) or the in-memory builder (which does
# support keep_fields/aux, via the "origin" index cpp_build_geom_graph
# returns mapping each output edge back to its input segment).
finalize_geom_graph <- function(x1, y1, x2, y2, cost, direction, snap_digits,
                                 attrs, keep_fields, output, out_prefix){
  if (output == "csv") {
    if (!is.null(keep_fields)) warning("keep_fields is ignored when output = 'csv' (that path never carries aux)")
    edges_path <- paste0(out_prefix, "_edges.csv")
    coords_path <- paste0(out_prefix, "_coords.csv")
    cpp_write_geom_graph_csv(x1, y1, x2, y2, cost, direction, snap_digits, edges_path, coords_path)
    return(invisible(c(edges_path, coords_path)))
  }

  result <- cpp_build_geom_graph(x1, y1, x2, y2, cost, direction, snap_digits)
  edges <- data.frame(from = result$from, to = result$to, dist = result$dist)
  coords <- data.frame(node_id = result$node_id, X = result$X, Y = result$Y)

  aux <- NULL
  if (!is.null(keep_fields)) {
    missing_fields <- setdiff(keep_fields, names(attrs))
    if (length(missing_fields) > 0) stop("keep_fields not found in attributes: ", paste(missing_fields, collapse = ", "))
    aux <- attrs[result$origin + 1L, keep_fields, drop = FALSE]
    rownames(aux) <- NULL
  }

  dispatch_graph_output(edges, coords, aux, output, out_prefix)
}

# Shared by every Module A/B reader: turn already-built edges/coords/aux
# tables into whichever of graph/data.frame/csv the caller asked for.
dispatch_graph_output <- function(edges, coords, aux, output, out_prefix){
  if (output == "data.frame") return(list(data = edges, coords = coords, aux = aux))

  if (output == "csv") {
    data.table::fwrite(edges, paste0(out_prefix, "_edges.csv"))
    data.table::fwrite(coords, paste0(out_prefix, "_coords.csv"))
    return(invisible(c(paste0(out_prefix, "_edges.csv"), paste0(out_prefix, "_coords.csv"))))
  }

  makegraph(edges, directed = TRUE, coords = coords, aux = aux)
}
