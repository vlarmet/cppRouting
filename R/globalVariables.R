# data.table's NSE (e.g. `:=`, `.N`, and column names used inside
# data.table expressions in graph_from_osmdata() and graph_from_sf())
# resolves column names as if they were variables. Declare them as
# global variables to satisfy R CMD check.
utils::globalVariables(c(
  ":=", ".N", "L1", "L2", "X", "Y",
  "coord_key", "feature", "id", "n_distinct", "osm_id", "part"
))
