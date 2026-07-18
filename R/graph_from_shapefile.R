#' Construct a graph from a shapefile of linestrings
#'
#' @param path Path to the \code{.shp} file (the companion \code{.dbf} must sit next to it, same base name).
#' @param cost_field Optional. Name of a numeric attribute column to use directly as edge cost.
#' @param speed_field Optional. Name of a numeric attribute column holding a speed. Used with \code{length_field}
#' (or the geometric length of each line, if \code{length_field} is \code{NULL}) to derive a cost as
#' \code{60 * length / speed}. Ignored if \code{cost_field} is provided.
#' @param length_field Optional. Name of a numeric attribute column holding each line's length. If \code{NULL}
#' and \code{speed_field} is used, length is computed from the line geometry itself.
#' @param direction_field Optional. Name of an attribute column describing the direction of travel allowed on
#' each line (e.g. a "sens" column). Combined with \code{direction_values}.
#' @param direction_values Optional. A list with up to three character vectors: \code{forward} (only from the
#' line's first point to its last), \code{backward} (only from last to first) and \code{twoway} (both directions).
#' Values of \code{direction_field} not listed anywhere default to \code{twoway}, with a warning.
#' @param directed logical. Used only when \code{direction_field} is \code{NULL}. If \code{TRUE}, every line is
#' kept as digitized (first point to last point) ; if \code{FALSE} (the default here), every line is made
#' bidirectional. Note this default intentionally differs from \link{makegraph}: without any direction
#' information, assuming every edge is traversable both ways is the safer default.
#' @param snap_digits integer. Number of decimal places used to round coordinates before deduplicating line
#' endpoints into nodes. Increase for very precise/projected coordinates, decrease if nearby endpoints coming
#' from slightly different digitizing are not being merged into the same node.
#' @param keep_fields Optional character vector of extra attribute columns to keep as \code{aux} (see
#' \link{makegraph}), e.g. an extra distance column alongside a time-based \code{cost_field}.
#' @param output One of \code{"graph"} (default, returns a ready-to-use \code{cppRouting_graph} via
#' \link{makegraph}), \code{"data.frame"} (returns the edges/coords tables without building the graph), or
#' \code{"csv"} (writes \code{<out_prefix>_edges.csv} and \code{<out_prefix>_coords.csv}).
#' @param out_prefix Path prefix used when \code{output = "csv"}.
#' @return Depending on \code{output}: a \code{cppRouting_graph}, a named list of two data.frames
#' (\code{data}, \code{coords}), or (invisibly) the paths of the written CSV files.
#' @details The shapefile is assumed to already be topologically clean, i.e. lines only meet at their
#' endpoints (the convention of reference road networks such as ROUTE500/BD TOPO). Junctions in the middle of
#' a line are not detected -- for that, see the OSM PBF import functions instead.
#' @export
graph_from_shapefile <- function(path,
                                  cost_field = NULL,
                                  speed_field = NULL,
                                  length_field = NULL,
                                  direction_field = NULL,
                                  direction_values = NULL,
                                  directed = FALSE,
                                  snap_digits = 7,
                                  keep_fields = NULL,
                                  output = c("graph", "data.frame", "csv"),
                                  out_prefix = NULL){

  output <- match.arg(output)
  if (output == "csv" && is.null(out_prefix)) stop("out_prefix must be provided when output = 'csv'")

  base <- sub("\\.shp$", "", path, ignore.case = TRUE)
  shp_path <- paste0(base, ".shp")
  dbf_path <- paste0(base, ".dbf")
  if (!file.exists(shp_path)) stop("Shapefile not found: ", shp_path)
  if (!file.exists(dbf_path)) stop("Companion .dbf not found: ", dbf_path)

  geom <- cpp_read_shp_geometry(shp_path)
  dbf <- cpp_read_dbf(dbf_path)
  if (nrow(dbf) != geom$n_records) {
    stop("Mismatch between number of records in .shp (", geom$n_records,
         ") and .dbf (", nrow(dbf), ")")
  }

  attrs <- dbf[geom$record + 1L, , drop = FALSE]

  prj_path <- paste0(base, ".prj")
  if (file.exists(prj_path) && is.null(length_field) && is.null(cost_field)) {
    prj_txt <- readLines(prj_path, warn = FALSE)
    if (grepl("GEOGCS", prj_txt[1], fixed = TRUE) && !grepl("PROJCS", prj_txt[1], fixed = TRUE)) {
      warning("Coordinates appear to be geographic (degrees), not projected. ",
              "Geometric length computed from these coordinates will not be a reliable planar distance; ",
              "consider reprojecting first or supplying length_field/cost_field directly.")
    }
  }

  if (!is.null(cost_field)) {
    if (!cost_field %in% names(attrs)) stop("cost_field '", cost_field, "' not found in .dbf attributes")
    cost <- as.numeric(attrs[[cost_field]])
  } else if (!is.null(speed_field)) {
    if (!speed_field %in% names(attrs)) stop("speed_field '", speed_field, "' not found in .dbf attributes")
    length_val <- if (!is.null(length_field)) {
      if (!length_field %in% names(attrs)) stop("length_field '", length_field, "' not found in .dbf attributes")
      as.numeric(attrs[[length_field]])
    } else geom$length
    cost <- 60 * length_val / as.numeric(attrs[[speed_field]])
  } else {
    cost <- geom$length
  }
  if (any(is.na(cost)) || any(cost < 0)) stop("Resolved cost contains NA or negative values")

  n <- length(geom$x1)
  direction_code <- rep(if (directed) 1L else 0L, n)
  if (!is.null(direction_field)) {
    if (!direction_field %in% names(attrs)) stop("direction_field '", direction_field, "' not found in .dbf attributes")
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
  }

  result <- cpp_build_geom_graph(geom$x1, geom$y1, geom$x2, geom$y2, cost, direction_code, snap_digits)

  edges <- data.frame(from = result$from, to = result$to, dist = result$dist)
  coords <- data.frame(node_id = result$node_id, X = result$X, Y = result$Y)

  aux <- NULL
  if (!is.null(keep_fields)) {
    missing_fields <- setdiff(keep_fields, names(attrs))
    if (length(missing_fields) > 0) stop("keep_fields not found in .dbf attributes: ", paste(missing_fields, collapse = ", "))
    aux <- attrs[result$origin + 1L, keep_fields, drop = FALSE]
    rownames(aux) <- NULL
  }

  if (output == "data.frame") return(list(data = edges, coords = coords, aux = aux))

  if (output == "csv") {
    data.table::fwrite(edges, paste0(out_prefix, "_edges.csv"))
    data.table::fwrite(coords, paste0(out_prefix, "_coords.csv"))
    return(invisible(c(paste0(out_prefix, "_edges.csv"), paste0(out_prefix, "_coords.csv"))))
  }

  makegraph(edges, directed = TRUE, coords = coords, aux = aux)
}
