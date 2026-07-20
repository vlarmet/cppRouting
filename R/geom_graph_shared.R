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

# A geographic (lat/lon degrees) coords table silently breaks the euclidean-
# distance heuristic the A*/NBA algorithms rely on (get_distance_pair()/
# get_path_pair() already document that they "require projected
# coordinates"): the heuristic ends up scaled in degrees while the edge cost
# is in meters/minutes/whatever, so A*/NBA either lose all speed benefit or,
# if a user tunes `constant` to compensate without truly reprojecting, can
# become inadmissible and return silently wrong (suboptimal) paths. Rather
# than let that trap ship silently, refuse outright unless the caller either
# doesn't want coordinates at all (coords = FALSE) or supplies epsg to
# reproject them into a projected CRS first.
stop_if_geographic_coords_unhandled <- function(is_geographic, want_coords, epsg){
  if (isTRUE(is_geographic) && want_coords && is.null(epsg)) {
    stop("Coordinates are geographic (lat/lon degrees), not projected. Euclidean distance between them does ",
         "not approximate the real edge cost, which silently breaks the A*/NBA algorithms' heuristic (see ",
         "get_distance_pair()/get_path_pair(): they require projected coordinates). Either set coords = FALSE ",
         "if you don't need node coordinates, or supply epsg = <a projected EPSG code> to reproject them ",
         "(requires the 'sf' package).")
  }
}

# Reprojects a node_id/X/Y coords data.frame from a geographic source CRS to
# a projected epsg via sf::st_transform() -- a no-op unless epsg is supplied
# AND the source coordinates are actually geographic (is_geographic), so sf
# is only ever required when reprojection is actually requested and needed.
# source_crs is anything sf::st_crs() accepts (an EPSG code, a WKT/PROJ4
# string, an existing crs object): shapefiles pass their raw .prj text
# (parsed directly by sf/PROJ, no EPSG lookup needed on our side), sf/osmdata
# sources pass sf::st_crs(x), OSM PBF sources hardcode 4326 (WGS84).
reproject_coords <- function(coords, is_geographic, source_crs, epsg){
  if (is.null(coords) || is.null(epsg) || !isTRUE(is_geographic)) return(coords)
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required to reproject coordinates (epsg argument). Install it, or leave epsg = NULL.")
  }
  pts <- sf::st_as_sf(coords, coords = c("X", "Y"), crs = source_crs)
  pts <- sf::st_transform(pts, epsg)
  xy <- sf::st_coordinates(pts)
  coords$X <- xy[, 1]
  coords$Y <- xy[, 2]
  coords
}

# Same, but operating on an already-written coords.csv in place -- used by
# the CSV-streaming writers (Module A/B), which write coordinates directly
# from C++ and never hold the full coords table in R otherwise.
reproject_coords_csv <- function(coords_path, is_geographic, source_crs, epsg){
  if (is.null(epsg) || !isTRUE(is_geographic)) return(invisible(NULL))
  coords <- data.table::fread(coords_path)
  coords <- reproject_coords(coords, is_geographic, source_crs, epsg)
  data.table::fwrite(coords, coords_path)
  invisible(NULL)
}

# Dispatches to either the CSV streaming writer (edges/coords written
# straight to disk as they're produced, aux never carried -- same
# restriction as the OSM CSV path) or the in-memory builder (which does
# support keep_fields/aux, via the "origin" index cpp_build_geom_graph
# returns mapping each output edge back to its input segment).
# want_coords: if FALSE, node coordinates are neither written (csv) nor
# attached to the returned graph/data.frame -- see cpp_extract's `coords`
# argument. Named want_coords here (not coords) because this function's own
# local variable holding the built coordinate table is already called
# `coords`.
finalize_geom_graph <- function(x1, y1, x2, y2, cost, direction, snap_digits,
                                 attrs, keep_fields, output, out_prefix, want_coords = TRUE,
                                 epsg = NULL, is_geographic = FALSE, source_crs = NULL){
  if (output == "csv") {
    if (!is.null(keep_fields)) warning("keep_fields is ignored when output = 'csv' (that path never carries aux)")
    edges_path <- paste0(out_prefix, "_edges.csv")
    coords_path <- paste0(out_prefix, "_coords.csv")
    cpp_write_geom_graph_csv(x1, y1, x2, y2, cost, direction, snap_digits, edges_path, coords_path, want_coords)
    if (want_coords) {
      reproject_coords_csv(coords_path, is_geographic, source_crs, epsg)
      return(invisible(c(edges_path, coords_path)))
    }
    return(invisible(edges_path))
  }

  result <- cpp_build_geom_graph(x1, y1, x2, y2, cost, direction, snap_digits)
  edges <- data.frame(from = result$from, to = result$to, dist = result$dist)
  coords <- data.frame(node_id = result$node_id, X = result$X, Y = result$Y)
  if (want_coords) coords <- reproject_coords(coords, is_geographic, source_crs, epsg)

  aux <- NULL
  if (!is.null(keep_fields)) {
    missing_fields <- setdiff(keep_fields, names(attrs))
    if (length(missing_fields) > 0) stop("keep_fields not found in attributes: ", paste(missing_fields, collapse = ", "))
    aux <- attrs[result$origin + 1L, keep_fields, drop = FALSE]
    rownames(aux) <- NULL
  }

  dispatch_graph_output(edges, coords, aux, output, out_prefix, want_coords)
}

# Shared by every Module A/B reader: turn already-built edges/coords/aux
# tables into whichever of graph/data.frame/csv the caller asked for.
# want_coords: if FALSE, coords is dropped (NULL in data.frame mode, not
# written in csv mode, not passed to makegraph() in graph mode).
dispatch_graph_output <- function(edges, coords, aux, output, out_prefix, want_coords = TRUE){
  if (!want_coords) coords <- NULL

  if (output == "data.frame") return(list(data = edges, coords = coords, aux = aux))

  if (output == "csv") {
    data.table::fwrite(edges, paste0(out_prefix, "_edges.csv"))
    if (!want_coords) return(invisible(paste0(out_prefix, "_edges.csv")))
    data.table::fwrite(coords, paste0(out_prefix, "_coords.csv"))
    return(invisible(c(paste0(out_prefix, "_edges.csv"), paste0(out_prefix, "_coords.csv"))))
  }

  makegraph(edges, directed = TRUE, coords = coords, aux = aux)
}
