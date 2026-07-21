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
#' (see \link{makegraph}). For \code{output = "csv"} (including the \code{max_memory} disk-backed
#' pipeline), this is written as an extra \code{<out_prefix>_aux.csv} instead, one row per edge in
#' \code{<out_prefix>_edges.csv} -- csv output always carries the same information as \code{data.frame}
#' output.
#' @param keep_length logical. If \code{TRUE}, keep each edge's haversine length (in km) as \code{aux},
#' regardless of the resolved time-based cost -- e.g. to get a time-based graph while still reporting
#' physical distance. Same \code{<out_prefix>_aux.csv} handling as \code{keep_clazz} for
#' \code{output = "csv"}.
#' @param output One of \code{"graph"} (default), \code{"data.frame"} or \code{"csv"}.
#' @param out_prefix Path prefix used when \code{output = "csv"}.
#' @param coords logical. If \code{TRUE} (default), extract each node's coordinates: written to
#' \code{<out_prefix>_coords.csv} (\code{output = "csv"}, including the \code{max_memory} disk-backed
#' pipeline), returned as a \code{data.frame(node_id, X, Y)} (\code{output = "data.frame"}), or passed to
#' \code{makegraph(coords = ...)} (\code{output = "graph"}). Set to \code{FALSE} to skip this -- for the
#' \code{max_memory} pipeline in particular, this also skips an entire external sort/merge pass over the
#' node coordinates, a real saving in disk I/O and time at continent scale, not just a smaller output. OSM
#' coordinates are always geographic (lat/lon), so \code{coords = TRUE} requires \code{epsg} to be set (see
#' below) -- otherwise this errors rather than silently returning coordinates that break the A*/NBA
#' algorithms.
#' @param epsg Optional. A projected EPSG code (or anything \code{sf::st_crs()} accepts) to reproject node
#' coordinates into, via \code{sf::st_transform()} (requires the \code{sf} package). Required whenever
#' \code{coords = TRUE} (OSM coordinates are always geographic). Not supported together with
#' \code{max_memory}: reprojecting would require reading the whole coords table back into R, defeating the
#' point of the disk-backed pipeline's bounded RAM -- reproject \verb{<out_prefix>_coords.csv} yourself
#' afterward instead (e.g. in chunks with \code{sf::st_transform()}).
#' @param max_memory Optional. A budget in megabytes for the disk-backed pipeline (see Details). When set,
#' \code{output} must be \code{"csv"}: the whole point is to never materialize the result (or the
#' intermediate node/way tables) beyond this budget, so returning a graph/data.frame to R would defeat it.
#' @param temp_dir Directory for the disk-backed pipeline's temporary files when \code{max_memory} is set.
#' Defaults to \code{tempdir()}.
#' @param progress logical. Print progress messages for each read pass.
#' @return Depending on \code{output}: a \code{cppRouting_graph}, a named list of two data.frames
#' (\code{data}, \code{coords}), or (invisibly) the paths of the written CSV files.
#' @details Unlike \link{graph_from_shapefile}/\link{graph_from_sf}, topology comes from OSM's own stable
#' node ids, not coordinate matching: a way can pass through a junction in its middle, and this is detected
#' (a node referenced by >= 2 kept ways, or a way's own endpoint) and split into separate edges accordingly.
#' A way referencing a node outside the file (a regional extract cut through a way) is dropped entirely, a
#' small, expected loss limited to the extract's edge. Turn restrictions, barriers and non-car profiles are
#' not supported yet.
#'
#' \code{output = "csv"} (without \code{max_memory}) streams edges/coordinates directly to disk as pass 3
#' produces them rather than building the full result in R memory first, keeping peak RAM bounded by passes
#' 1/2's in-RAM structures alone -- suitable for country-scale extracts (a few GB peak RAM even for all of
#' France).
#'
#' \code{max_memory} goes further: instead of holding the referenced-node table and kept ways in RAM at all
#' (passes 1/2 above), every intermediate table is spilled to disk via external (chunked) sort and
#' sequential merge-joins, keeping peak RAM close to the stated budget regardless of the input file's total
#' size -- meant for files that don't fit in RAM even after the passes-1/2-only optimization above (e.g. a
#' full-planet extract, or a country-scale extract on a memory-constrained machine). This is slower (more
#' disk I/O, several sort/merge passes) and requires nodes to appear in strictly increasing id order in the
#' file, which is standard for osmium/osmconvert/Geofabrik extracts but not guaranteed by the PBF format
#' itself -- violated files raise a clear error rather than silently producing a wrong graph.
#' @noRd
graph_from_osmpbf <- function(path,
                               profile = osm_car_profile(),
                               use_maxspeed = TRUE,
                               assume_oneway_motorway = TRUE,
                               keep_clazz = FALSE,
                               keep_length = FALSE,
                               output = c("graph", "data.frame", "csv"),
                               out_prefix = NULL,
                               coords = TRUE,
                               epsg = NULL,
                               max_memory = NULL,
                               temp_dir = tempdir(),
                               progress = TRUE){

  output <- match.arg(output)
  if (output == "csv" && is.null(out_prefix)) stop("out_prefix must be provided when output = 'csv'")
  if (!is.null(max_memory)) {
    if (output != "csv") {
      stop("max_memory requires output = 'csv': the disk-backed pipeline never materializes the full result ",
           "in R memory, by design -- read the CSV files back yourself afterward if you know they'll fit.")
    }
    if (!is.numeric(max_memory) || length(max_memory) != 1 || is.na(max_memory) || max_memory <= 0) {
      stop("max_memory must be a single positive number (megabytes), not: ", paste(deparse(max_memory), collapse = " "))
    }
    if (!is.null(epsg) && coords) {
      stop("epsg is not supported together with max_memory: reprojecting would require reading the full ",
           "coords table back into R, defeating the whole point of the disk-backed pipeline's bounded RAM. ",
           "Either drop max_memory, set coords = FALSE, or reproject <out_prefix>_coords.csv yourself ",
           "afterward (e.g. in chunks with sf::st_transform()).")
    }
  }
  stop_if_geographic_coords_unhandled(TRUE, coords, epsg)
  if (!file.exists(path)) stop("File not found: ", path)
  if (!all(c("highway", "clazz", "speed") %in% names(profile))) {
    stop("profile must have columns 'highway', 'clazz', 'speed'")
  }

  if (!is.null(max_memory)) {
    edges_path <- paste0(out_prefix, "_edges.csv")
    coords_path <- paste0(out_prefix, "_coords.csv")
    aux_path <- paste0(out_prefix, "_aux.csv")
    cpp_build_osm_graph_disk_csv(path,
                                  as.character(profile$highway),
                                  as.integer(profile$clazz),
                                  as.numeric(profile$speed),
                                  use_maxspeed, assume_oneway_motorway,
                                  as.numeric(max_memory), temp_dir,
                                  edges_path, coords_path, aux_path,
                                  coords, keep_clazz, keep_length, progress)
    paths <- edges_path
    if (coords) paths <- c(paths, coords_path)
    if (keep_clazz || keep_length) paths <- c(paths, aux_path)
    return(invisible(paths))
  }

  if (output == "csv") {
    # Streams edges/coords/aux straight to disk during the same junction-
    # splitting pass instead of building the full result in R memory first --
    # keeps peak RAM bounded by pass 1/2's structures alone, for
    # continent-scale files.
    edges_path <- paste0(out_prefix, "_edges.csv")
    coords_path <- paste0(out_prefix, "_coords.csv")
    aux_path <- paste0(out_prefix, "_aux.csv")
    cpp_write_osm_graph_csv(path,
                             as.character(profile$highway),
                             as.integer(profile$clazz),
                             as.numeric(profile$speed),
                             use_maxspeed, assume_oneway_motorway,
                             edges_path, coords_path, aux_path,
                             coords, keep_clazz, keep_length, progress)
    if (coords) reproject_coords_csv(coords_path, TRUE, 4326L, epsg)
    paths <- edges_path
    if (coords) paths <- c(paths, coords_path)
    if (keep_clazz || keep_length) paths <- c(paths, aux_path)
    return(invisible(paths))
  }

  result <- cpp_build_osm_graph(path,
                                 as.character(profile$highway),
                                 as.integer(profile$clazz),
                                 as.numeric(profile$speed),
                                 use_maxspeed, assume_oneway_motorway, progress)

  edges <- data.frame(from = result$from, to = result$to, dist = result$dist)
  coords_df <- data.frame(node_id = result$node_id, X = result$lon, Y = result$lat)
  if (coords) coords_df <- reproject_coords(coords_df, TRUE, 4326L, epsg)
  aux <- NULL
  if (keep_clazz) aux <- data.frame(clazz = result$clazz)
  if (keep_length) {
    if (is.null(aux)) aux <- data.frame(length = result$length) else aux$length <- result$length
  }

  dispatch_graph_output(edges, coords_df, aux, output, out_prefix, want_coords = coords)
}
