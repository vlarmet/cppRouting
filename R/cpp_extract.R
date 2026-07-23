#' Extract a routable graph from a geographic source
#'
#' Single entry point to build a \link{makegraph}-ready graph from a geographic linestring source. The
#' source is dispatched on automatically from \code{x}'s type (file extension, or object class). Exactly
#' five forms of \code{x} are recognized -- see the \code{x} parameter below for the precise contract of
#' each, since passing the wrong shape of object for a given source (e.g. only part of an \code{osmdata}
#' result) fails loudly rather than silently producing a wrong graph.
#'
#' Arguments are grouped by which family of source they apply to: arguments that work the same way
#' regardless of source stay as their own top-level parameter; source-specific ones are grouped into the
#' \code{osm} or \code{non_osm} list (see those two parameters below for the full list of what each one
#' accepts) -- so a shapefile/\code{sf} call never needs to mention \code{profile}/\code{use_maxspeed}, and
#' a \code{.pbf}/\code{osmdata}/\code{osmextract} call never needs to mention \code{cost_field}/
#' \code{direction_field}.
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
#'   cached \code{.pbf}} with \code{osmextract::oe_get(..., download_only = TRUE, skip_vectortranslate =
#'   TRUE)} (\strong{both} flags are required: \code{download_only} alone still runs the GDAL conversion and
#'   returns a \code{.gpkg} path instead) and pass that \code{.pbf} path like any other \code{.pbf} file.
#'   Do not let osmextract convert it to \code{sf}/GeoPackage via GDAL first: that conversion only keeps
#'   tagged nodes by default, which is not enough to reconstruct topology.
#' }
#' GeoPackage is not supported and has no plans to be: convert it to a shapefile first (e.g. via
#' \code{sf::st_write(sf::st_read("file.gpkg"), "file.shp")}, or with any GIS tool) and use that instead.
#'
#' @param output One of \code{"graph"} (default, returns a ready-to-use \code{cppRouting_graph}),
#' \code{"data.frame"}, or \code{"csv"}. Common to every source.
#' @param out_prefix Path prefix used when \code{output = "csv"}. Common to every source.
#' @param coords logical. If \code{TRUE} (default), extract each node's coordinates: written to
#' \code{<out_prefix>_coords.csv} (\code{output = "csv"}, including the \code{.pbf} \code{max_memory}
#' disk-backed pipeline), returned as a \code{data.frame(node_id, X, Y)} (\code{output = "data.frame"}), or
#' passed to \code{makegraph(coords = ...)} (\code{output = "graph"}). Set to \code{FALSE} to skip this --
#' e.g. if you already have node coordinates elsewhere and only need the graph topology. For the
#' \code{max_memory} pipeline this also skips an entire external sort/merge pass, a real saving in disk I/O
#' and time at continent scale, not just a smaller output. If the source coordinates are geographic
#' (lat/lon -- always true for \code{.pbf}/\code{osmdata}) and \code{coords = TRUE}, \code{epsg} must be
#' supplied (see below), otherwise this errors rather than silently returning coordinates that break the
#' A*/NBA algorithms (see \link{get_distance_pair}). Common to every source.
#' @param epsg Optional. A projected EPSG code (or anything \code{sf::st_crs()} accepts) to reproject node
#' coordinates into, via \code{sf::st_transform()} (requires the \code{sf} package). Only takes effect if
#' \code{coords = TRUE} and the source CRS is geographic; already-projected coordinates (possible for
#' shapefile/sf sources) are left as-is regardless. Not supported together with \code{osm$max_memory} on a
#' \code{.pbf} source: reprojecting would require reading the whole coords table back into R, defeating the
#' point of the disk-backed pipeline's bounded RAM -- reproject \verb{<out_prefix>_coords.csv} yourself
#' afterward instead. Common to every source.
#' @param keep_length logical. If \code{TRUE}, keep each edge's length (in km) as \code{aux}, regardless of
#' how the main \code{cost} was actually resolved -- e.g. to get a time-based graph while still reporting
#' physical distance. Shapefile/sf sources: the geometric length (unreliable for a geographic CRS, same
#' caveat as \code{non_osm$length_field} below -- reproject or use \code{non_osm$keep_fields} with your own
#' length column instead in that case). OSM sources: the haversine length. \code{output = "csv"} carries
#' the same information as \code{output = "data.frame"}: \code{aux} is written to \verb{<out_prefix>_aux.csv}.
#' Common to every source.
#' @param snap_digits integer. Shapefile/sf sources: decimal places used to round coordinates before
#' deduplicating line endpoints into nodes. \code{osmdata} sources: decimal places used to match line
#' vertex coordinates against \code{$osm_points} coordinates. Unused for \code{.pbf} (topology comes from
#' OSM node ids directly, no coordinate matching needed). Defaults (\code{NULL}) to \code{7} for geographic
#' (lat/lon) coordinates or \code{4} for projected (meter-scale) ones -- geographic coordinates need many
#' decimal places to distinguish nearby points, but that precision is far too tight for projected
#' coordinates, where it starts picking up floating-point noise between independently-digitized but
#' logically-identical points rather than real geometric differences. The one argument that doesn't sort
#' cleanly into \code{osm}/\code{non_osm}: it applies to shapefile/sf sources and to \code{osmdata}, but not
#' to \code{.pbf}.
#'
#' @param osm A list of options for OSM sources (\code{.pbf}/\code{osmdata}/\code{osmextract}) only, ignored
#' for shapefile/\code{sf}. Any subset may be supplied -- omitted elements keep their default below; an
#' unrecognized element name raises an error immediately (rather than silently keeping its default, which is
#' what a typo'd list element would otherwise do). Recognized elements:
#' \itemize{
#'   \item \code{mode}: one of \code{"car"} (default), \code{"foot"}, \code{"bicycle"} or \code{"hgv"} --
#'   which access-tag/oneway resolution rules to apply (see \link{graph_from_osmpbf}'s \code{@details} for
#'   the exact tags checked per mode). If \code{profile} is left \code{NULL}, also selects the matching
#'   default profile (\link{osm_car_profile}/\link{osm_foot_profile}/\link{osm_bike_profile}/
#'   \link{osm_hgv_profile}).
#'   \item \code{profile}: a data.frame with columns \code{highway}, \code{clazz}, \code{speed} (see
#'   \code{osm_car_profile} and friends): which \code{highway} tag values to keep, and their class id /
#'   default speed (km/h) when not overridden by \code{maxspeed}. Defaults (\code{NULL}) to the profile
#'   matching \code{mode}; pass your own to override independently of \code{mode}.
#'   \item \code{use_maxspeed}: logical, default \code{TRUE}. If \code{TRUE}, a way's \code{maxspeed} tag
#'   overrides the profile's default speed when present and parseable (plain km/h, or \code{"<n> mph"}). No
#'   effect for \code{mode = "foot"}/\code{"bicycle"}.
#'   \item \code{assume_oneway_motorway}: logical, default \code{TRUE}. If \code{TRUE}, \code{motorway}/
#'   \code{motorway_link} ways are treated as one-way unless explicitly tagged \code{oneway=no}. Irrelevant
#'   for \code{mode = "foot"} (oneway is always ignored for that mode).
#'   \item \code{keep_clazz}: logical, default \code{FALSE}. If \code{TRUE}, keep each edge's resolved
#'   highway class id as \code{aux}.
#'   \item \code{max_memory}: \code{.pbf} only (ignored for \code{osmdata}). A budget in megabytes for a
#'   disk-backed pipeline that spills every intermediate table (not just the final result) to disk, keeping
#'   peak RAM close to this budget regardless of the input file's size -- for files too large to fit in RAM
#'   even with the default pipeline (e.g. a full-planet extract). Requires \code{output = "csv"}. See
#'   \code{graph_from_osmpbf} for details and the (rare but possible) requirement that nodes appear in
#'   increasing id order in the file.
#'   \item \code{temp_dir}: \code{.pbf} only. Directory for temporary files when \code{max_memory} is set.
#'   Defaults to \code{tempdir()}.
#'   \item \code{progress}: logical, default \code{TRUE}. Print progress messages while reading/matching.
#' }
#' @param non_osm A list of options for shapefile/\code{sf} sources only, ignored for OSM sources. Any
#' subset may be supplied -- omitted elements keep their default below; an unrecognized element name raises
#' an error immediately (same reasoning as \code{osm} above). Recognized elements:
#' \itemize{
#'   \item \code{cost_field}: optional, name of a numeric attribute column to use directly as edge cost.
#'   \item \code{speed_field}: optional, name of a numeric attribute column holding a speed in km/h. Used
#'   with \code{length_field} (in km) to derive a cost in minutes as \code{60 * length / speed}. Ignored if
#'   \code{cost_field} is provided.
#'   \item \code{length_field}: optional, name of a numeric attribute column holding each line's length, in
#'   km (to match the \code{60 * length / speed} formula above). If \code{NULL} and \code{speed_field} is
#'   used, length is instead computed from the line geometry itself and converted from the CRS's native unit
#'   (assumed to be meters, true for any typical projected CRS) to km.
#'   \item \code{direction_field}: optional, name of an attribute column describing the direction of travel
#'   allowed on each line (e.g. a "sens" column). Combined with \code{direction_values}.
#'   \item \code{direction_values}: optional, a list with up to three character vectors: \code{forward}
#'   (only from the line's first point to its last), \code{backward} (only from last to first) and
#'   \code{twoway} (both directions). Values of \code{direction_field} not listed anywhere default to
#'   \code{twoway}, with a warning.
#'   \item \code{directed}: logical, default \code{FALSE}. Used only when \code{direction_field} is
#'   \code{NULL}. If \code{TRUE}, every line is kept as digitized (first point to last point); if
#'   \code{FALSE}, every line is made bidirectional.
#'   \item \code{keep_fields}: optional character vector of extra attribute columns to keep as \code{aux}
#'   (see \link{makegraph}).
#' }
#'
#' @return Depending on \code{output}: a \code{cppRouting_graph}, a named list of two data.frames
#' (\code{data}, \code{coords}), or (invisibly) the paths of the written CSV files. For
#' \code{output = "data.frame"}, \code{from}/\code{to}/\code{node_id} are \code{character} regardless of
#' source: a node id is a label, not a quantity, matching what \code{makegraph()} itself already guarantees
#' for \code{dict$ref}. For \code{output = "graph"}, that same real id is carried through \code{dict$ref};
#' \code{data$from}/\code{data$to} stay \code{makegraph()}'s own compact integer index (its adjacency list is
#' always built on that, not on the original id).
#' @details Shapefile/sf sources are assumed to already be topologically clean, i.e. lines only meet at
#' their endpoints. OSM sources
#' (\code{.pbf} or \code{osmdata}) instead detect junctions in the middle of a way from OSM's own node ids
#' (a node referenced by 2 or more kept ways, or a way's own endpoint); for \code{.pbf}, a way referencing a
#' node outside the file (a regional extract cut through it) is dropped entirely, and for \code{osmdata}, a
#' way with a vertex that can't be matched against \code{$osm_points} is dropped, with a warning. A plain
#' \code{sf} object that looks OSM-derived (it has an \code{osm_id} column) but isn't the full
#' \code{osmdata} result is refused with an explanatory error, rather than silently treated as
#' topologically clean like a plain \code{sf} source and producing a wrong graph.
#' @export
cpp_extract <- function(x,
                         output = c("graph", "data.frame", "csv"),
                         out_prefix = NULL,
                         coords = TRUE,
                         epsg = NULL,
                         keep_length = FALSE,
                         snap_digits = NULL,
                         osm = list(),
                         non_osm = list()){

  output <- match.arg(output)

  osm_defaults <- list(mode = "car", profile = NULL, use_maxspeed = TRUE,
                        assume_oneway_motorway = TRUE, keep_clazz = FALSE,
                        max_memory = NULL, temp_dir = tempdir(), progress = TRUE)
  non_osm_defaults <- list(cost_field = NULL, speed_field = NULL, length_field = NULL,
                            direction_field = NULL, direction_values = NULL,
                            directed = FALSE, keep_fields = NULL)

  osm <- merge_arg_list(osm, osm_defaults, "osm")
  non_osm <- merge_arg_list(non_osm, non_osm_defaults, "non_osm")

  args <- c(non_osm, list(snap_digits = snap_digits, keep_length = keep_length,
                           output = output, out_prefix = out_prefix, coords = coords, epsg = epsg))

  osm_args <- list(mode = osm$mode, profile = osm$profile, use_maxspeed = osm$use_maxspeed,
                    assume_oneway_motorway = osm$assume_oneway_motorway,
                    keep_clazz = osm$keep_clazz, keep_length = keep_length, output = output,
                    out_prefix = out_prefix, coords = coords, epsg = epsg, progress = osm$progress)

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
           "instead, don't let it convert to sf/gpkg via GDAL: get the raw cached .pbf path with ",
           "osmextract::oe_get(..., download_only = TRUE, skip_vectortranslate = TRUE) (both flags are ",
           "required) and pass that .pbf path to cpp_extract() directly -- treating this sf object as a ",
           "plain topologically-clean network (like ",
           "graph_from_sf() does) would silently produce a wrong graph, since OSM ways can have ",
           "junctions in the middle of a line.")
    }
    return(do.call(graph_from_sf, c(list(x = x), args)))
  }

  if (is.character(x) && length(x) == 1) {
    ext <- tolower(tools::file_ext(x))
    if (ext == "shp") return(do.call(graph_from_shapefile, c(list(path = x), args)))
    if (ext == "gpkg") {
      stop("GeoPackage is not supported. Convert it to a shapefile first, e.g. ",
           "sf::st_write(sf::st_read('file.gpkg'), 'file.shp'), and pass that .shp path instead.")
    }
    if (ext == "pbf") {
      return(do.call(graph_from_osmpbf, c(list(path = x, max_memory = osm$max_memory, temp_dir = osm$temp_dir), osm_args)))
    }
    if (ext == "osm") stop("Plain OSM XML (.osm) is not supported, only the binary .osm.pbf format")
    stop("Unrecognized file extension '.", ext, "'. Supported: .shp, .pbf")
  }

  stop("cpp_extract() does not know how to handle an object of class: ", paste(class(x), collapse = "/"))
}
