#' Construct a graph from an OSM PBF file
#'
#' Internal implementation for the \code{.pbf}/\code{.osm.pbf} source of \link{cpp_extract}; not exported
#' directly.
#' @param path Path to the \code{.osm.pbf} file.
#' @param profile A data.frame with columns \code{highway}, \code{clazz}, \code{speed} (see
#' \link{osm_car_profile} for the default): which \code{highway} tag values to keep, and their class id /
#' default speed (km/h) when not overridden by \code{maxspeed}.
#' @param use_maxspeed logical. If \code{TRUE} (default), a way's \code{maxspeed} tag overrides the
#' profile's default speed when present and parseable (plain km/h, or \code{"<n> mph"}).
#' @param assume_oneway_motorway logical. If \code{TRUE} (default), \code{motorway}/\code{motorway_link}
#' ways are treated as one-way (forward, first node to last) unless an explicit \code{oneway=no} tag says
#' otherwise. Standard OSM \code{oneway} tag values (\code{yes}/\code{-1}/\code{no}) always take precedence.
#' @param keep_clazz logical. If \code{TRUE}, keep each edge's resolved highway class id as \code{aux}
#' (see \link{makegraph}).
#' @param output One of \code{"graph"} (default), \code{"data.frame"} or \code{"csv"}.
#' @param out_prefix Path prefix used when \code{output = "csv"}.
#' @param progress logical. Print progress messages for each of the 3 read passes.
#' @return Depending on \code{output}: a \code{cppRouting_graph}, a named list of two data.frames
#' (\code{data}, \code{coords}), or (invisibly) the paths of the written CSV files.
#' @details Unlike \link{graph_from_shapefile}/\link{graph_from_sf}, topology comes from OSM's own stable
#' node ids, not coordinate matching: a way can pass through a junction in its middle, and this is detected
#' (a node referenced by >= 2 kept ways, or a way's own endpoint) and split into separate edges accordingly.
#' A way referencing a node outside the file (a regional extract cut through a way) is dropped entirely, a
#' small, expected loss limited to the extract's edge. Turn restrictions, barriers and non-car profiles are
#' not supported yet.
#' @noRd
graph_from_osmpbf <- function(path,
                               profile = osm_car_profile(),
                               use_maxspeed = TRUE,
                               assume_oneway_motorway = TRUE,
                               keep_clazz = FALSE,
                               output = c("graph", "data.frame", "csv"),
                               out_prefix = NULL,
                               progress = TRUE){

  output <- match.arg(output)
  if (output == "csv" && is.null(out_prefix)) stop("out_prefix must be provided when output = 'csv'")
  if (!file.exists(path)) stop("File not found: ", path)
  if (!all(c("highway", "clazz", "speed") %in% names(profile))) {
    stop("profile must have columns 'highway', 'clazz', 'speed'")
  }

  result <- cpp_build_osm_graph(path,
                                 as.character(profile$highway),
                                 as.integer(profile$clazz),
                                 as.numeric(profile$speed),
                                 use_maxspeed, assume_oneway_motorway, progress)

  edges <- data.frame(from = result$from, to = result$to, dist = result$dist)
  coords <- data.frame(node_id = result$node_id, X = result$lon, Y = result$lat)
  aux <- if (keep_clazz) data.frame(clazz = result$clazz) else NULL

  dispatch_graph_output(edges, coords, aux, output, out_prefix)
}
