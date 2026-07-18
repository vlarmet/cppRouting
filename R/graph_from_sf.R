#' Construct a graph from an in-memory sf object of linestrings
#'
#' Internal implementation for the \code{sf} object source of \link{cpp_extract}; not exported directly.
#' @param x An \code{sf} object with \code{LINESTRING}/\code{MULTILINESTRING} geometries (e.g. one already
#' read via \code{sf::st_read()}, or coming from \code{osmdata}/\code{osmextract}).
#' @inheritParams graph_from_shapefile
#' @return Depending on \code{output}: a \code{cppRouting_graph}, a named list of two data.frames
#' (\code{data}, \code{coords}), or (invisibly) the paths of the written CSV files.
#' @details Same options and topology assumption as \link{graph_from_shapefile}: \code{x} is assumed to
#' already be topologically clean, i.e. lines only meet at their endpoints. \code{MULTILINESTRING} features
#' are split into one edge per part. Requires the \code{sf} package.
#' @noRd
graph_from_sf <- function(x,
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
  if (!requireNamespace("sf", quietly = TRUE)) stop("Package 'sf' is required for graph_from_sf()")
  if (!inherits(x, "sf")) stop("x must be an sf object")
  if (output == "csv" && is.null(out_prefix)) stop("out_prefix must be provided when output = 'csv'")

  geom_type <- as.character(sf::st_geometry_type(x, by_geometry = TRUE))
  if (!all(geom_type %in% c("LINESTRING", "MULTILINESTRING"))) {
    stop("graph_from_sf() only supports LINESTRING/MULTILINESTRING geometries, found: ",
         paste(unique(geom_type), collapse = ", "))
  }

  # st_coordinates() already returns a part id (L1) and, for MULTILINESTRING,
  # a separate feature id (L2) -- for a homogeneous input this needs no
  # per-feature branching at all. st_cast() would give the same shape
  # uniformly but scales badly (roughly quadratically) on real-sized inputs,
  # so mixed LINESTRING/MULTILINESTRING inputs are instead split into their
  # two homogeneous subsets and processed separately.
  extract_xy <- function(x_sub, orig_index){
    xy <- data.table::as.data.table(sf::st_coordinates(x_sub))
    if ("L2" %in% names(xy)) {
      xy[, `:=`(feature = orig_index[L2], part = L1)]
    } else {
      xy[, `:=`(feature = orig_index[L1], part = 1L)]
    }
    xy
  }

  is_multi <- geom_type == "MULTILINESTRING"
  xy_list <- list()
  if (any(!is_multi)) xy_list[[length(xy_list) + 1]] <- extract_xy(x[!is_multi, ], which(!is_multi))
  if (any(is_multi)) xy_list[[length(xy_list) + 1]] <- extract_xy(x[is_multi, ], which(is_multi))
  xy <- data.table::rbindlist(xy_list, fill = TRUE)

  parts <- xy[, list(x1 = X[1L], y1 = Y[1L], x2 = X[.N], y2 = Y[.N],
                      length = sum(sqrt(diff(X)^2 + diff(Y)^2))),
              by = list(feature, part)]
  data.table::setorder(parts, feature, part)

  attrs <- sf::st_drop_geometry(x)[parts$feature, , drop = FALSE]

  is_geographic <- isTRUE(suppressMessages(suppressWarnings(sf::st_is_longlat(x))))
  warn_if_geographic_without_length(is_geographic, length_field, cost_field)

  cost <- resolve_geom_cost(attrs, parts$length, cost_field, speed_field, length_field)
  direction_code <- resolve_geom_direction(attrs, direction_field, direction_values, directed, nrow(parts))

  result <- cpp_build_geom_graph(parts$x1, parts$y1, parts$x2, parts$y2, cost, direction_code, snap_digits)

  finalize_geom_graph(result, attrs, keep_fields, output, out_prefix)
}
