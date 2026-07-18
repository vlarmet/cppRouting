.onLoad <- function (libname, pkgname) {

  Sys.setenv ("RCPP_PARALLEL_BACKEND" = "tinythread")
  invisible ()
}

# Required for data.table's [.data.table NSE (column names resolved as if
# they were variables, e.g. in graph_from_sf()'s by-group aggregation) to
# work correctly when called from inside this package's namespace -- without
# it, cedta() falls back to a conservative mode and column lookups like
# `X`/`Y`/`L1`/`L2` fail with "object not found".
.datatable.aware <- TRUE
