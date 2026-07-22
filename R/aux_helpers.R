# Internal helpers shared by get_distance_pair(), get_distance_matrix() (and,
# eventually, other query functions) to support aggregate_aux with either a
# single legacy aux vector or a multi-column aux data.frame declared at
# makegraph() time.

# Turn (aux attribute, aggregate_aux argument) into a numeric matrix ready for
# the C++ layer, plus enough metadata to shape the result the same way the
# request was expressed.
resolve_aggregate_aux <- function(aux, aggregate_aux){
  if (is.null(aux)) stop("network don't have additional weight")
  aux_is_df <- is.data.frame(aux)

  if (isTRUE(aggregate_aux)) {
    cols <- if (aux_is_df) colnames(aux) else "aux"
  } else {
    cols <- as.character(aggregate_aux)
    if (length(cols) == 0) stop("aggregate_aux must be TRUE or a character vector of aux column names")
  }

  if (aux_is_df) {
    missing <- setdiff(cols, colnames(aux))
    if (length(missing) > 0) stop(sprintf("aux column(s) not found: %s", paste(missing, collapse = ", ")))
    mat <- as.matrix(aux[, cols, drop = FALSE])
  } else {
    unknown <- setdiff(cols, "aux")
    if (length(unknown) > 0) stop(sprintf("aux column(s) not found: %s", paste(unknown, collapse = ", ")))
    mat <- matrix(as.double(aux), ncol = 1, dimnames = list(NULL, "aux"))
  }
  storage.mode(mat) <- "double"

  list(matrix = mat, names = cols, single_legacy = isTRUE(aggregate_aux) && !aux_is_df)
}

# Reshape a (n_pairs x n_cols) numeric matrix of aggregated values back into a
# plain vector for the single-column legacy case (aggregate_aux = TRUE on a
# graph with one aux column, matching the historical return shape), or a
# data.frame with named columns otherwise.
shape_aggregate_result <- function(value_matrix, col_names, single_legacy){
  if (single_legacy) return(as.numeric(value_matrix[, 1]))
  out <- as.data.frame(value_matrix)
  colnames(out) <- col_names
  out
}
