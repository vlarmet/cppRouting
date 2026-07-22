# Internal helpers shared by graph_from_shapefile(), graph_from_geopackage()
# and graph_from_sf() (Module A: "geometry endpoint graph builder", see the
# graph-import plan). Each format-specific reader only has to produce a flat
# per-segment table (x1,y1,x2,y2,length,attrs) -- everything from cost/
# direction resolution down to the final graph/data.frame/csv output is
# common and lives here.

# Merges a user-supplied list of overrides against a named list of defaults,
# used by cpp_extract()'s `osm`/`non_osm` grouped-argument lists. A flat
# function argument gets typo protection for free from R's own "unused
# argument" error; a list element does not (out$us_maxspeed <- FALSE would
# just silently add an ignored element, leaving use_maxspeed at its default)
# -- so unrecognized names are rejected explicitly here instead.
merge_arg_list <- function(user_list, defaults, arg_name){
  if (is.null(user_list)) return(defaults)
  if (!is.list(user_list)) stop("'", arg_name, "' must be a list")
  unknown <- setdiff(names(user_list), names(defaults))
  if (length(unknown) > 0) {
    stop("'", arg_name, "' has unrecognized element(s): ", paste(unknown, collapse = ", "),
         ". Valid elements are: ", paste(names(defaults), collapse = ", "))
  }
  utils::modifyList(defaults, user_list)
}

# snap_digits = NULL resolves to a CRS-appropriate default rather than one
# flat number: geographic (lat/lon degree) coordinates need many decimal
# places to distinguish nearby points (7 decimals ~= 1cm), but that same
# precision is far too tight for projected (meter-scale) coordinates like
# Lambert-93 -- real-world validation against IGN's ROUTE500 NOEUD_ROUTIER
# layer showed exact-match rates collapsing from 99.9% at snap_digits = 4 to
# 0% at snap_digits = 6, purely from floating-point/production noise between
# independently-stored coordinates of the "same" point (observed mismatches
# were ~5 micrometers apart -- nowhere near a real geometric difference).
resolve_snap_digits <- function(snap_digits, is_geographic){
  if (!is.null(snap_digits)) return(snap_digits)
  if (isTRUE(is_geographic)) 7L else 4L
}

resolve_geom_cost <- function(attrs, geom_length, cost_field, speed_field, length_field, is_geographic){
  if (!is.null(cost_field)) {
    if (!cost_field %in% names(attrs)) stop("cost_field '", cost_field, "' not found in attributes")
    cost <- as.numeric(attrs[[cost_field]])
  } else if (!is.null(speed_field)) {
    if (!speed_field %in% names(attrs)) stop("speed_field '", speed_field, "' not found in attributes")
    length_val <- if (!is.null(length_field)) {
      if (!length_field %in% names(attrs)) stop("length_field '", length_field, "' not found in attributes")
      as.numeric(attrs[[length_field]])
    } else if (!isTRUE(is_geographic)) {
      # geom_length is computed directly from the (projected, meter-scale)
      # coordinates -- convert to km so that 60*length/speed gives minutes
      # for a speed in km/h, same convention as an explicit length_field
      # (e.g. ROUTE500's LONGUEUR, itself in km). Without this, relying on
      # speed_field without an explicit length_field silently produced a
      # cost 1000x too large. Not applied for geographic coordinates
      # (degrees): warn_if_geographic_without_length() already flags that
      # combination as unreliable regardless.
      geom_length / 1000
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
  # node_id round-trips through fread() above as a plain double (the C++
  # writer already appended ".0" for any id past the int32 range -- see
  # avoid_integer64_ambiguity()) -- but fwrite() below would print a
  # whole-number double back out as plain digits, silently re-introducing
  # the exact bit64::integer64 ambiguity this rewrite pass would otherwise
  # undo. Reapply it before writing the final file.
  coords$node_id <- avoid_integer64_ambiguity(as.character(coords$node_id))
  data.table::fwrite(coords, coords_path)
  invisible(NULL)
}

# Dispatches to either the CSV streaming writer (edges/coords written
# straight to disk as they're produced -- the fast path, used whenever no
# aux is requested) or the in-memory builder (needed to get the "origin"
# index cpp_build_geom_graph returns, mapping each output edge back to its
# input segment, for keep_fields/keep_length). CSV output always carries the
# same information as data.frame output -- an extra "<out_prefix>_aux.csv"
# is written whenever keep_fields/keep_length produced one, exactly mirroring
# data.frame mode's separate `aux` element (see dispatch_graph_output).
# want_coords: if FALSE, node coordinates are neither written (csv) nor
# attached to the returned graph/data.frame -- see cpp_extract's `coords`
# argument. Named want_coords here (not coords) because this function's own
# local variable holding the built coordinate table is already called
# `coords`.
finalize_geom_graph <- function(x1, y1, x2, y2, cost, direction, snap_digits,
                                 attrs, keep_fields, output, out_prefix, want_coords = TRUE,
                                 epsg = NULL, is_geographic = FALSE, source_crs = NULL,
                                 geom_length = NULL, keep_length = FALSE){
  needs_origin <- !is.null(keep_fields) || keep_length
  if (output == "csv" && !needs_origin) {
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
  if (keep_length) {
    # Independent of how cost itself was resolved (cost_field/speed_field/
    # fallback): always the geometric length, converted from the CRS's
    # native unit (assumed meters) to km -- except for geographic sources,
    # where it's left in degrees (unreliable as a physical distance either
    # way; see warn_if_geographic_without_length()).
    length_km <- geom_length[result$origin + 1L]
    if (!isTRUE(is_geographic)) length_km <- length_km / 1000
    if (is.null(aux)) aux <- data.frame(length = length_km) else aux$length <- length_km
  }

  dispatch_graph_output(edges, coords, aux, output, out_prefix, want_coords)
}

# from/to/node_id are character (node ids are labels, not quantities -- see
# cpp_extract()'s docs), but a large real OSM id (see graph_from_osmdata(),
# the one Module B path whose ids reach dispatch_graph_output() rather than
# a dedicated C++ CSV writer) is still a string of plain digits, and
# fread()'s type inference does not consult the source type: it sees plain
# digits and, once they exceed 2^31, promotes the column to
# bit64::integer64 regardless of quoting -- which without the bit64 package
# installed (a Suggests-only dependency of neither package) holds outright
# garbage (an int64 bit pattern reinterpreted as a plain double), not merely
# a display quirk. Appending ".0" to the digit string (not a numeric
# round-trip -- the original text is kept, just extended) makes the column
# unambiguously a double on read, exact for any id up to 2^53 -- comfortably
# past any real OSM id. Only applied when actually needed (values beyond
# int32 range), so Module A's small dense ids -- never at risk -- keep their
# existing, unaffected plain-digit CSV formatting.
avoid_integer64_ambiguity <- function(x){
  if (!is.character(x)) return(x)
  vals <- suppressWarnings(as.numeric(x))
  needs_fix <- !is.na(vals) & abs(vals) > .Machine$integer.max
  x[needs_fix] <- paste0(x[needs_fix], ".0")
  x
}

# Shared by every Module A/B reader: turn already-built edges/coords/aux
# tables into whichever of graph/data.frame/csv the caller asked for.
# want_coords: if FALSE, coords is dropped (NULL in data.frame mode, not
# written in csv mode, not passed to makegraph() in graph mode). csv mode
# writes the exact same tables data.frame mode would return, one file per
# table (aux only when non-NULL) -- there is no reason for csv and
# data.frame to disagree on what the result is, only on how it's returned.
dispatch_graph_output <- function(edges, coords, aux, output, out_prefix, want_coords = TRUE){
  if (!want_coords) coords <- NULL

  if (output == "data.frame") return(list(data = edges, coords = coords, aux = aux))

  if (output == "csv") {
    edges_out <- edges
    edges_out$from <- avoid_integer64_ambiguity(edges_out$from)
    edges_out$to <- avoid_integer64_ambiguity(edges_out$to)
    data.table::fwrite(edges_out, paste0(out_prefix, "_edges.csv"))
    paths <- paste0(out_prefix, "_edges.csv")
    if (want_coords) {
      coords_out <- coords
      coords_out$node_id <- avoid_integer64_ambiguity(coords_out$node_id)
      data.table::fwrite(coords_out, paste0(out_prefix, "_coords.csv"))
      paths <- c(paths, paste0(out_prefix, "_coords.csv"))
    }
    if (!is.null(aux)) {
      data.table::fwrite(aux, paste0(out_prefix, "_aux.csv"))
      paths <- c(paths, paste0(out_prefix, "_aux.csv"))
    }
    return(invisible(paths))
  }

  # cpp_extract()'s aux is a named data.frame of extra per-edge attributes
  # (keep_fields/keep_clazz/keep_length columns) -- passed straight through,
  # since makegraph()'s aux natively accepts a named data.frame with one or
  # more columns (see R/graphs.R), keeping each column's name and carrying
  # all of them into the graph's attrib$aux at once.
  makegraph(edges, directed = TRUE, coords = coords, aux = aux)
}
