#' Extract a routable graph from a geographic source
#'
#' Single entry point to build a \link{makegraph}-ready graph from a geographic linestring source. The
#' source is dispatched on automatically from \code{x}'s type (file extension, or object class). Exactly
#' five forms of \code{x} are recognized -- see the \code{x} parameter below for the precise contract of
#' each, since passing the wrong shape of object for a given source (e.g. only part of an \code{osmdata}
#' result) fails loudly rather than silently producing a wrong graph.
#'
#' @param x One of:
#' \itemize{
#'   \item a path (character) to a \strong{\code{.shp}} shapefile -- the companion \code{.dbf}/\code{.prj}
#'   must sit next to it, same base name;
#'   \item an in-memory \strong{\code{sf}} object (\code{LINESTRING}/\code{MULTILINESTRING}), e.g. from
#'   \code{sf::st_read()} on a GeoPackage or any other OGR-readable format;
#'   \item a path (character) to an \strong{\code{.osm.pbf}} file;
#'   \item the \strong{full result of \code{osmdata::osmdata_sf()}} (class \code{"osmdata"}) -- not just its
#'   \code{$osm_lines} element. \code{$osm_points} must also be present: it holds the \code{osm_id} of every
#'   node osmdata downloaded (including untagged vertices in the middle of a way), which is what lets
#'   \code{cpp_extract} recover real OSM topology and detect junctions;
#'   \item for \strong{\code{osmextract}}: not a distinct input type -- get the path to the \strong{raw
#'   cached \code{.pbf}} (e.g. \code{osmextract::oe_get(..., download_only = TRUE)}) and pass that path like
#'   any other \code{.pbf} file. Do not let osmextract convert it to \code{sf}/GeoPackage via GDAL first:
#'   that conversion only keeps tagged nodes by default, which is not enough to reconstruct topology.
#' }
#' GeoPackage paths are not implemented yet (deferred).
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
#' bidirectional.
#' @param snap_digits integer. Shapefile/sf sources: decimal places used to round coordinates before
#' deduplicating line endpoints into nodes. \code{osmdata} sources: decimal places used to match line
#' vertex coordinates against \code{$osm_points} coordinates. Unused for \code{.pbf} (topology comes from
#' OSM node ids directly, no coordinate matching needed).
#' @param keep_fields Optional character vector of extra attribute columns to keep as \code{aux} (see
#' \link{makegraph}).
#' @param output One of \code{"graph"} (default, returns a ready-to-use \code{cppRouting_graph}),
#' \code{"data.frame"}, or \code{"csv"}.
#' @param out_prefix Path prefix used when \code{output = "csv"}.
#' @param profile OSM sources (\code{.pbf}/\code{osmdata}) only. A data.frame with columns \code{highway},
#' \code{clazz}, \code{speed} describing which \code{highway} tag values to keep and their class id /
#' default speed (km/h). Defaults to \link{osm_car_profile}.
#' @param use_maxspeed OSM sources only. If \code{TRUE} (default), a way's \code{maxspeed} tag overrides the
#' profile's default speed when present and parseable.
#' @param assume_oneway_motorway OSM sources only. If \code{TRUE} (default), \code{motorway}/
#' \code{motorway_link} ways are treated as one-way unless explicitly tagged \code{oneway=no}.
#' @param keep_clazz OSM sources only. If \code{TRUE}, keep each edge's resolved highway class id as
#' \code{aux}.
#' @param progress OSM sources only. Print progress messages while reading/matching.
#' @return Depending on \code{output}: a \code{cppRouting_graph}, a named list of two data.frames
#' (\code{data}, \code{coords}), or (invisibly) the paths of the written CSV files.
#' @details Shapefile/sf sources are assumed to already be topologically clean, i.e. lines only meet at
#' their endpoints (the convention of reference road networks such as ROUTE500/BD TOPO). OSM sources
#' (\code{.pbf} or \code{osmdata}) instead detect junctions in the middle of a way from OSM's own node ids
#' (a node referenced by 2 or more kept ways, or a way's own endpoint); for \code{.pbf}, a way referencing a
#' node outside the file (a regional extract cut through it) is dropped entirely, and for \code{osmdata}, a
#' way with a vertex that can't be matched against \code{$osm_points} is dropped, with a warning. A plain
#' \code{sf} object that looks OSM-derived (it has an \code{osm_id} column) but isn't the full
#' \code{osmdata} result is refused with an explanatory error, rather than silently treated as
#' topologically clean like a plain \code{sf} source and producing a wrong graph.
#' @export
cpp_extract <- function(x,
                         cost_field = NULL,
                         speed_field = NULL,
                         length_field = NULL,
                         direction_field = NULL,
                         direction_values = NULL,
                         directed = FALSE,
                         snap_digits = 7,
                         keep_fields = NULL,
                         output = c("graph", "data.frame", "csv"),
                         out_prefix = NULL,
                         profile = osm_car_profile(),
                         use_maxspeed = TRUE,
                         assume_oneway_motorway = TRUE,
                         keep_clazz = FALSE,
                         progress = TRUE){

  output <- match.arg(output)
  args <- list(cost_field = cost_field, speed_field = speed_field, length_field = length_field,
               direction_field = direction_field, direction_values = direction_values,
               directed = directed, snap_digits = snap_digits, keep_fields = keep_fields,
               output = output, out_prefix = out_prefix)

  osm_args <- list(profile = profile, use_maxspeed = use_maxspeed,
                    assume_oneway_motorway = assume_oneway_motorway,
                    keep_clazz = keep_clazz, output = output,
                    out_prefix = out_prefix, progress = progress)

  # Class is the reliable signal for a full osmdata result (class "osmdata",
  # a list -- not itself an "sf" object): route it to the osmdata adapter,
  # which needs BOTH $osm_lines and $osm_points (checked inside
  # graph_from_osmdata) to reconstruct real OSM topology.
  if (inherits(x, "osmdata")) {
    return(do.call(graph_from_osmdata, c(list(x = x, snap_digits = snap_digits), osm_args)))
  }

  if (inherits(x, "sf")) {
    if ("osm_id" %in% names(x)) {
      stop("This sf object has an 'osm_id' column but isn't the full osmdata result (class 'osmdata'). ",
           "Topology reconstruction needs $osm_points alongside $osm_lines -- pass the whole result of ",
           "osmdata::osmdata_sf() instead of just its $osm_lines element. If this came from osmextract ",
           "instead, don't let it convert to sf/gpkg via GDAL: get the raw cached .pbf path (e.g. ",
           "osmextract::oe_get(..., download_only = TRUE)) and pass that .pbf path to cpp_extract() ",
           "directly -- treating this sf object as a plain topologically-clean network (like ",
           "graph_from_sf() does) would silently produce a wrong graph, since OSM ways can have ",
           "junctions in the middle of a line.")
    }
    return(do.call(graph_from_sf, c(list(x = x), args)))
  }

  if (is.character(x) && length(x) == 1) {
    ext <- tolower(tools::file_ext(x))
    if (ext == "shp") return(do.call(graph_from_shapefile, c(list(path = x), args)))
    if (ext == "gpkg") stop("GeoPackage sources are not implemented yet")
    if (ext == "pbf") return(do.call(graph_from_osmpbf, c(list(path = x), osm_args)))
    if (ext == "osm") stop("Plain OSM XML (.osm) is not supported, only the binary .osm.pbf format")
    stop("Unrecognized file extension '.", ext, "'. Supported: .shp, .pbf")
  }

  stop("cpp_extract() does not know how to handle an object of class: ", paste(class(x), collapse = "/"))
}
