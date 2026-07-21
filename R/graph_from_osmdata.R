#' Construct a graph from an osmdata result
#'
#' Internal implementation for the \code{osmdata} source of \link{cpp_extract}; not exported directly.
#' @param x The full result of \code{osmdata::osmdata_sf()} (class \code{"osmdata"}), with both
#' \code{$osm_lines} and \code{$osm_points} present. Passing only \code{x$osm_lines} does not work: OSM
#' topology (which nodes are junctions) can only be recovered by matching every vertex of every line
#' against \code{$osm_points}, which holds the \code{osm_id} of every node osmdata downloaded, including
#' untagged vertices in the middle of a way. \code{$osm_lines} alone has no such per-vertex id.
#' @inheritParams graph_from_osmpbf
#' @param snap_digits integer. Decimal places used when matching line vertex coordinates against
#' \code{$osm_points} coordinates. Defaults (\code{NULL}) to \code{7} for geographic (lat/lon) coordinates --
#' the normal case, since \code{osmdata} results are essentially always WGS84 -- or \code{4} for the rare
#' case of an already-projected \code{osmdata} object; see \link{graph_from_shapefile}'s \code{snap_digits}
#' for why the two need different precision.
#' @return Depending on \code{output}: a \code{cppRouting_graph}, a named list of two data.frames
#' (\code{data}, \code{coords}), or (invisibly) the paths of the written CSV files.
#' @details Every non-\code{osm_id} column of \code{x$osm_lines} is treated as an OSM tag (its column name
#' as the key, its value for that row as the value) when resolving the profile/access/oneway/maxspeed
#' rules -- exactly matching how \code{osmdata_sf()} flattens tags into columns. A line whose vertex cannot
#' be matched to any \code{$osm_points} coordinate (beyond \code{snap_digits} tolerance) is dropped, with a
#' warning.
#' @noRd
graph_from_osmdata <- function(x,
                                profile = osm_car_profile(),
                                use_maxspeed = TRUE,
                                assume_oneway_motorway = TRUE,
                                keep_clazz = FALSE,
                                keep_length = FALSE,
                                snap_digits = NULL,
                                output = c("graph", "data.frame", "csv"),
                                out_prefix = NULL,
                                coords = TRUE,
                                epsg = NULL,
                                progress = TRUE){

  output <- match.arg(output)
  if (!inherits(x, "osmdata")) stop("x must be the full result of osmdata::osmdata_sf() (class 'osmdata')")
  if (is.null(x$osm_lines) || is.null(x$osm_points)) {
    stop("x must contain both $osm_lines and $osm_points. Topology reconstruction needs $osm_points ",
         "(the osm_id of every node, including untagged vertices in the middle of a way) -- passing only ",
         "$osm_lines discards exactly the information needed to detect junctions.")
  }
  if (!requireNamespace("sf", quietly = TRUE)) stop("Package 'sf' is required for graph_from_osmdata()")
  if (output == "csv" && is.null(out_prefix)) stop("out_prefix must be provided when output = 'csv'")

  lines <- x$osm_lines
  points <- x$osm_points

  is_geographic <- isTRUE(suppressMessages(suppressWarnings(sf::st_is_longlat(points))))
  stop_if_geographic_coords_unhandled(is_geographic, coords, epsg)
  snap_digits <- resolve_snap_digits(snap_digits, is_geographic)

  geom_type <- as.character(sf::st_geometry_type(lines, by_geometry = TRUE))
  if (!all(geom_type %in% c("LINESTRING", "MULTILINESTRING"))) {
    stop("graph_from_osmdata() only supports LINESTRING/MULTILINESTRING geometries in x$osm_lines, found: ",
         paste(unique(geom_type), collapse = ", "))
  }

  key <- function(x, y) sprintf("%.*f_%.*f", snap_digits, round(x, snap_digits), snap_digits, round(y, snap_digits))

  pt_xy <- sf::st_coordinates(points)
  pt_key <- key(pt_xy[, 1], pt_xy[, 2])

  # A coordinate key shared by two DIFFERENT osm_id values (not the same
  # node repeated) means the lookup below can only keep one of them --
  # whichever vertex matching happens to land on that key silently gets
  # the wrong node id, with no error. This only bites when snap_digits is
  # too coarse for how densely packed the actual points are (a real risk
  # in dense urban areas), so it's worth surfacing rather than staying quiet.
  collisions <- data.table::data.table(coord_key = pt_key, id = as.character(points$osm_id))[
    , list(n_distinct = length(unique(id))), by = coord_key][n_distinct > 1]
  if (nrow(collisions) > 0) {
    warning(nrow(collisions), " coordinate location(s) in $osm_points are shared by multiple distinct osm_id ",
            "values at snap_digits = ", snap_digits, " precision -- only one id is kept at each such location, ",
            "which can silently misattribute any way vertex that lands there to the wrong node. Increase ",
            "snap_digits if this count seems too high for the data's actual density.")
  }

  pt_id <- stats::setNames(as.character(points$osm_id), pt_key)

  extract_refs <- function(lines_sub, orig_index){
    xy <- data.table::as.data.table(sf::st_coordinates(lines_sub))
    if ("L2" %in% names(xy)) xy[, `:=`(feature = orig_index[L2], part = L1)]
    else xy[, `:=`(feature = orig_index[L1], part = 1L)]
    xy[, osm_id := pt_id[key(X, Y)]]
    xy
  }

  is_multi <- geom_type == "MULTILINESTRING"
  xy_list <- list()
  if (any(!is_multi)) xy_list[[length(xy_list) + 1]] <- extract_refs(lines[!is_multi, ], which(!is_multi))
  if (any(is_multi)) xy_list[[length(xy_list) + 1]] <- extract_refs(lines[is_multi, ], which(is_multi))
  xy <- data.table::rbindlist(xy_list, fill = TRUE)

  data.table::setorder(xy, feature, part)
  ways <- xy[, list(refs = list(osm_id)), by = list(feature, part)]

  incomplete <- vapply(ways$refs, anyNA, logical(1))
  if (any(incomplete)) {
    warning(sum(incomplete), " way part(s) had a vertex that could not be matched to an osm_id in ",
            "$osm_points (beyond snap_digits tolerance, or the point was filtered out of the query) -- ",
            "dropped.")
    ways <- ways[!incomplete]
  }
  if (nrow(ways) == 0) stop("No way could be matched against $osm_points -- nothing to build a graph from")

  attrs <- sf::st_drop_geometry(lines)[ways$feature, , drop = FALSE]
  tag_cols <- setdiff(names(attrs), "osm_id")

  n <- nrow(attrs)
  way_tag_keys <- vector("list", n)
  way_tag_vals <- vector("list", n)
  for (i in seq_len(n)) {
    vals <- as.character(unlist(attrs[i, tag_cols], use.names = FALSE))
    keep <- !is.na(vals) & vals != ""
    way_tag_keys[[i]] <- tag_cols[keep]
    way_tag_vals[[i]] <- vals[keep]
  }

  if (progress) message("[cpp_extract] matched ", n, " way part(s) against ", nrow(points), " osm_points; resolving profile and building edges...")

  result <- cpp_build_osm_graph_from_data(lapply(ways$refs, as.numeric), way_tag_keys, way_tag_vals,
                                           as.numeric(points$osm_id), pt_xy[, 2], pt_xy[, 1],
                                           as.character(profile$highway), as.integer(profile$clazz), as.numeric(profile$speed),
                                           use_maxspeed, assume_oneway_motorway)

  edges <- data.frame(from = result$from, to = result$to, dist = result$dist)
  coords_df <- data.frame(node_id = result$node_id, X = result$lon, Y = result$lat)
  if (coords) coords_df <- reproject_coords(coords_df, is_geographic, sf::st_crs(points), epsg)
  aux <- NULL
  if (keep_clazz) aux <- data.frame(clazz = result$clazz)
  if (keep_length) {
    if (is.null(aux)) aux <- data.frame(length = result$length) else aux$length <- result$length
  }

  dispatch_graph_output(edges, coords_df, aux, output, out_prefix, want_coords = coords)
}
