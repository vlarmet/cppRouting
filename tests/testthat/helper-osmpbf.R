# Minimal hand-rolled OSM PBF *encoder*, used only to build tiny test
# fixtures (the package itself only ever needs to decode). Mirrors exactly
# the subset of fileformat.proto/osmformat.proto that src/io_osmpbf.cpp
# decodes: BlobHeader/Blob (raw, uncompressed -- no need to exercise zlib
# here, that's exercised implicitly by every real .osm.pbf test against
# Bourgogne/France), StringTable, DenseNodes, Way.

pbf_varint <- function(v){
  v <- as.numeric(v)
  out <- integer(0)
  repeat {
    b <- bitwAnd(v %% 128, 0x7F)
    v <- v %/% 128
    if (v > 0) b <- bitwOr(b, 0x80)
    out <- c(out, b)
    if (v <= 0) break
  }
  as.raw(out)
}

pbf_zigzag <- function(v) ifelse(v >= 0, 2 * v, -2 * v - 1)

pbf_tag <- function(field_no, wire_type) pbf_varint(bitwOr(bitwShiftL(field_no, 3), wire_type))

pbf_bytes_field <- function(field_no, raw_bytes) c(pbf_tag(field_no, 2), pbf_varint(length(raw_bytes)), raw_bytes)

pbf_varint_field <- function(field_no, v) c(pbf_tag(field_no, 0), pbf_varint(v))

pbf_packed_svarint_field <- function(field_no, values){
  body <- unlist(lapply(pbf_zigzag(values), pbf_varint))
  pbf_bytes_field(field_no, as.raw(body))
}

pbf_packed_varint_field <- function(field_no, values){
  body <- unlist(lapply(values, pbf_varint))
  pbf_bytes_field(field_no, as.raw(body))
}

pbf_string_field <- function(field_no, s) pbf_bytes_field(field_no, charToRaw(s))

# nodes: data.frame(id, lat, lon). ways: list of list(refs=int vector, tags=named list of key=value strings).
write_test_osmpbf <- function(path, nodes, ways){
  strings <- unique(unlist(c(list(""), lapply(ways, function(w) c(names(w$tags), unlist(w$tags))))))
  str_index <- setNames(seq_along(strings) - 1L, strings)

  stringtable <- unlist(lapply(strings, pbf_string_field, field_no = 1))

  dense_ids <- c(diff(c(0, nodes$id)))
  dense_lats <- diff(c(0, round(nodes$lat * 1e9 / 100)))
  dense_lons <- diff(c(0, round(nodes$lon * 1e9 / 100)))
  dense_nodes <- c(pbf_packed_svarint_field(1, dense_ids),
                    pbf_packed_svarint_field(8, dense_lats),
                    pbf_packed_svarint_field(9, dense_lons))

  way_bytes <- lapply(seq_along(ways), function(i){
    w <- ways[[i]]
    keys <- if (length(w$tags)) unname(str_index[names(w$tags)]) else integer(0)
    vals <- if (length(w$tags)) unname(str_index[unlist(w$tags)]) else integer(0)
    ref_deltas <- diff(c(0, w$refs))
    body <- c(pbf_varint_field(1, i),
              if (length(keys)) pbf_packed_varint_field(2, keys) else raw(0),
              if (length(vals)) pbf_packed_varint_field(3, vals) else raw(0),
              pbf_packed_svarint_field(8, ref_deltas))
    body
  })

  group <- c(pbf_bytes_field(2, as.raw(dense_nodes)),
             unlist(lapply(way_bytes, function(b) pbf_bytes_field(3, as.raw(b)))))

  primitive_block <- c(pbf_bytes_field(1, as.raw(stringtable)),
                        pbf_bytes_field(2, as.raw(group)))

  blob <- pbf_bytes_field(1, as.raw(primitive_block)) # raw (uncompressed) payload

  blob_header <- c(pbf_string_field(1, "OSMData"), pbf_varint_field(3, length(blob)))

  header_len <- length(blob_header)
  len_be <- as.raw(c(bitwAnd(bitwShiftR(header_len, 24), 0xFF), bitwAnd(bitwShiftR(header_len, 16), 0xFF),
                      bitwAnd(bitwShiftR(header_len, 8), 0xFF), bitwAnd(header_len, 0xFF)))

  con <- file(path, "wb")
  on.exit(close(con))
  writeBin(as.raw(c(len_be, blob_header, blob)), con)
}

# Small synthetic network shared by test-osmpbf.R and test-osm_disk_pipeline.R
# (must live here, not in a test- file, so it's available regardless of
# alphabetical test file ordering): three nodes on one way (A-B-C, id 1-2-3)
# plus a second way starting at the *middle* node of the first one (B-D, id
# 2-4). Node 2 is therefore a true mid-way junction -- exactly the case
# graph_from_shapefile() cannot handle (endpoint-only topology) but the whole
# point of the OSM PBF path is to detect it and split way A-B-C into two edges.
make_test_pbf <- function(){
  path <- tempfile(fileext = ".pbf")
  nodes <- data.frame(id = 1:4, lat = c(47.00, 47.01, 47.02, 47.01), lon = c(5.00, 5.00, 5.00, 5.01))
  ways <- list(
    list(refs = c(1, 2, 3), tags = list(highway = "residential")),
    list(refs = c(2, 4), tags = list(highway = "residential", oneway = "yes"))
  )
  write_test_osmpbf(path, nodes, ways)
  path
}
