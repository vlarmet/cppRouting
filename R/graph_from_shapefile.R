#' Construct a graph from a shapefile of linestrings
#'
#' Internal implementation for the \code{.shp} source of \link{cpp_extract}; not exported directly.
#' @param path Path to the \code{.shp} file (the companion \code{.dbf} must sit next to it, same base name).
#' @param cost_field Optional. Name of a numeric attribute column to use directly as edge cost.
#' @param speed_field Optional. Name of a numeric attribute column holding a speed in km/h. Used with
#' \code{length_field} (in km) to derive a cost in minutes as \code{60 * length / speed}. Ignored if
#' \code{cost_field} is provided.
#' @param length_field Optional. Name of a numeric attribute column holding each line's length, in km (to
#' match the \code{60 * length / speed} formula above). If \code{NULL} and \code{speed_field} is used,
#' length is instead computed from the line geometry itself and converted from the CRS's native unit
#' (assumed to be meters, true for any typical projected CRS) to km.
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
#' endpoints into nodes. Defaults (\code{NULL}) to \code{7} for a geographic (lat/lon) CRS or \code{4} for a
#' projected one, based on the \code{.prj} file -- geographic coordinates need many decimal places to
#' distinguish nearby points, but that same precision is far too tight for projected (meter-scale)
#' coordinates, where it starts picking up floating-point noise between independently-digitized but
#' logically-identical points rather than real geometric differences. Decrease further if nearby endpoints
#' coming from slightly different digitizing are still not being merged into the same node.
#' @param keep_fields Optional character vector of extra attribute columns to keep as \code{aux} (see
#' \link{makegraph}), e.g. an extra distance column alongside a time-based \code{cost_field}.
#' @param keep_length logical. If \code{TRUE}, keep each edge's geometric length (in km) as \code{aux},
#' regardless of how \code{cost} was actually resolved (\code{cost_field}/\code{speed_field}/fallback) --
#' e.g. to get a time-based graph while still reporting physical distance. For a geographic CRS, this is
#' the raw degree-based length (as unreliable as a physical distance as the geometric-length cost fallback
#' -- see \code{length_field} above); reproject (\code{epsg}) or use \code{keep_fields} with your own length
#' column instead in that case. \code{output = "csv"} carries the same information as
#' \code{output = "data.frame"}: \code{aux} is written to \verb{<out_prefix>_aux.csv}.
#' @param output One of \code{"graph"} (default, returns a ready-to-use \code{cppRouting_graph} via
#' \link{makegraph}), \code{"data.frame"} (returns the edges/coords tables without building the graph), or
#' \code{"csv"} (writes \code{<out_prefix>_edges.csv} and \code{<out_prefix>_coords.csv}).
#' @param out_prefix Path prefix used when \code{output = "csv"}.
#' @param coords logical. If \code{TRUE} (default), extract each node's coordinates: written to
#' \code{<out_prefix>_coords.csv} (\code{output = "csv"}), returned as a \code{data.frame(node_id, X, Y)}
#' (\code{output = "data.frame"}), or passed to \code{makegraph(coords = ...)} (\code{output = "graph"}). Set
#' to \code{FALSE} to skip this -- e.g. if you already have node coordinates elsewhere and only need the
#' graph topology. If the shapefile's \code{.prj} declares a geographic (lat/lon) CRS and \code{coords} is
#' \code{TRUE}, \code{epsg} must be supplied (see below) -- otherwise this errors rather than silently
#' returning coordinates that break the A*/NBA algorithms.
#' @param epsg Optional. A projected EPSG code (or anything \code{sf::st_crs()} accepts) to reproject node
#' coordinates into, via \code{sf::st_transform()} (requires the \code{sf} package). Only takes effect if
#' \code{coords = TRUE} and the source CRS is geographic; otherwise ignored (already-projected coordinates
#' are left as-is). \code{A*}/\code{NBA} (see \link{get_distance_pair}) need coordinates in the same
#' (projected) unit as the edge cost -- geographic (lat/lon) coordinates silently break their heuristic.
#' @return Depending on \code{output}: a \code{cppRouting_graph}, a named list of two data.frames
#' (\code{data}, \code{coords}), or (invisibly) the paths of the written CSV files.
#' @details The shapefile is assumed to already be topologically clean, i.e. lines only meet at their
#' endpoints (the convention of reference road networks such as ROUTE500/BD TOPO). Junctions in the middle of
#' a line are not detected -- for that, see the OSM PBF import functions instead.
#' @noRd
graph_from_shapefile <- function(path,
                                  cost_field = NULL,
                                  speed_field = NULL,
                                  length_field = NULL,
                                  direction_field = NULL,
                                  direction_values = NULL,
                                  directed = FALSE,
                                  snap_digits = NULL,
                                  keep_fields = NULL,
                                  keep_length = FALSE,
                                  output = c("graph", "data.frame", "csv"),
                                  out_prefix = NULL,
                                  coords = TRUE,
                                  epsg = NULL){

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
  is_geographic <- FALSE
  source_crs <- NULL
  if (file.exists(prj_path)) {
    prj_txt <- readLines(prj_path, warn = FALSE)
    is_geographic <- grepl("GEOGCS", prj_txt[1], fixed = TRUE) && !grepl("PROJCS", prj_txt[1], fixed = TRUE)
    source_crs <- paste(prj_txt, collapse = " ")
  }
  warn_if_geographic_without_length(is_geographic, length_field, cost_field)
  stop_if_geographic_coords_unhandled(is_geographic, coords, epsg)
  snap_digits <- resolve_snap_digits(snap_digits, is_geographic)

  cost <- resolve_geom_cost(attrs, geom$length, cost_field, speed_field, length_field, is_geographic)
  direction_code <- resolve_geom_direction(attrs, direction_field, direction_values, directed, length(geom$x1))

  finalize_geom_graph(geom$x1, geom$y1, geom$x2, geom$y2, cost, direction_code, snap_digits,
                       attrs, keep_fields, output, out_prefix, want_coords = coords,
                       epsg = epsg, is_geographic = is_geographic, source_crs = source_crs,
                       geom_length = geom$length, keep_length = keep_length)
}
