#' Extract a routable graph from a geographic source
#'
#' Single entry point to build a \link{makegraph}-ready graph from a geographic linestring source: a
#' shapefile path, an in-memory \code{sf} object, or (as more source types are added) a GeoPackage path,
#' an OSM PBF file, or an \code{osmdata}/\code{osmextract} result. The source is dispatched on automatically
#' from \code{x}'s type (file extension, or object class).
#'
#' @param x A path to a source file (currently \code{.shp}), or an in-memory \code{sf} object.
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
#' @param snap_digits integer. Number of decimal places used to round coordinates before deduplicating line
#' endpoints into nodes.
#' @param keep_fields Optional character vector of extra attribute columns to keep as \code{aux} (see
#' \link{makegraph}).
#' @param output One of \code{"graph"} (default, returns a ready-to-use \code{cppRouting_graph}),
#' \code{"data.frame"}, or \code{"csv"}.
#' @param out_prefix Path prefix used when \code{output = "csv"}.
#' @param profile OSM \code{.pbf} sources only. A data.frame with columns \code{highway}, \code{clazz},
#' \code{speed} describing which \code{highway} tag values to keep and their class id / default speed
#' (km/h). Defaults to \link{osm_car_profile}.
#' @param use_maxspeed OSM \code{.pbf} sources only. If \code{TRUE} (default), a way's \code{maxspeed} tag
#' overrides the profile's default speed when present and parseable.
#' @param assume_oneway_motorway OSM \code{.pbf} sources only. If \code{TRUE} (default),
#' \code{motorway}/\code{motorway_link} ways are treated as one-way unless explicitly tagged \code{oneway=no}.
#' @param keep_clazz OSM \code{.pbf} sources only. If \code{TRUE}, keep each edge's resolved highway class
#' id as \code{aux}.
#' @param progress OSM \code{.pbf} sources only. Print progress messages for each of the 3 read passes.
#' @return Depending on \code{output}: a \code{cppRouting_graph}, a named list of two data.frames
#' (\code{data}, \code{coords}), or (invisibly) the paths of the written CSV files.
#' @details Sources of type shapefile/geopackage/sf are assumed to already be topologically clean, i.e.
#' lines only meet at their endpoints (the convention of reference road networks such as ROUTE500/BD TOPO).
#' OSM \code{.pbf} sources instead detect junctions in the middle of a way from OSM's own node ids (a node
#' referenced by 2 or more kept ways, or a way's own endpoint), and a way referencing a node outside the
#' file (a regional extract cut through it) is dropped entirely. \code{cpp_extract} refuses (rather than
#' silently mishandling) an object it
#' recognizes as OSM-derived but not a raw \code{.pbf} file -- either the full result of
#' \code{osmdata::osmdata_sf()} (class \code{"osmdata"}), or a plain \code{sf} object with an \code{osm_id}
#' column (as returned by \code{osmdata$osm_lines} or \code{osmextract}).
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

  osm_source_error <- function(reason){
    stop("This looks like an OSM-derived object (", reason, "). Unlike a plain reference-network sf ",
         "object, OSM ways can have junctions in the middle of a line, not just at its endpoints -- ",
         "treating it as topologically clean (like graph_from_sf() does) would silently produce a wrong ",
         "graph. OSM-aware handling is not implemented yet; use a raw .pbf source once OSM support lands, ",
         "or pre-split ways at shared nodes yourself before calling cpp_extract().")
  }

  # Class is the reliable signal for a full osmdata result (class "osmdata",
  # a list -- not itself an "sf" object). osmextract and a manually-extracted
  # osmdata$osm_lines are both plain "sf" objects with no distinguishing
  # class, so those can only be caught by the 'osm_id' column both give their
  # feature id column, matching OSM's own way/node id.
  if (inherits(x, "osmdata")) osm_source_error("class 'osmdata'")

  if (inherits(x, "sf")) {
    if ("osm_id" %in% names(x)) osm_source_error("it has an 'osm_id' column")
    return(do.call(graph_from_sf, c(list(x = x), args)))
  }

  if (is.character(x) && length(x) == 1) {
    ext <- tolower(tools::file_ext(x))
    if (ext == "shp") return(do.call(graph_from_shapefile, c(list(path = x), args)))
    if (ext == "gpkg") stop("GeoPackage sources are not implemented yet")
    if (ext == "pbf") {
      return(graph_from_osmpbf(x, profile = profile, use_maxspeed = use_maxspeed,
                                assume_oneway_motorway = assume_oneway_motorway,
                                keep_clazz = keep_clazz, output = output,
                                out_prefix = out_prefix, progress = progress))
    }
    if (ext == "osm") stop("Plain OSM XML (.osm) is not supported, only the binary .osm.pbf format")
    stop("Unrecognized file extension '.", ext, "'. Supported: .shp, .pbf")
  }

  stop("cpp_extract() does not know how to handle an object of class: ", paste(class(x), collapse = "/"))
}
