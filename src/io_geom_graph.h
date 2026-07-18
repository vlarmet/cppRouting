#ifndef IO_GEOM_GRAPH_H
#define IO_GEOM_GRAPH_H

#include <Rcpp.h>

// Module A: build a compact from/to/dist edge table (plus a node_id/X/Y
// coordinate table) from a flat list of segment endpoints, by deduplicating
// endpoints that fall on the same coordinate (rounded to snap_digits decimal
// places). This is the generalization of the endpoint-matching approach used
// for ROUTE500-like networks (data_readme/data_preparation.R): topology is
// carried entirely by coincident endpoint coordinates, never by a shared
// node id (those don't exist for shapefile/geopackage/sf sources).
//
// direction[i] convention: 0 = twoway, 1 = forward only (x1,y1 -> x2,y2),
// 2 = backward only (x2,y2 -> x1,y1). Self-loops (same node on both ends
// after snapping) are dropped, matching data_preparation.R's own filter.
// The returned "origin" vector maps each output edge row back to its input
// segment index (0-based), so callers can replicate per-segment attributes
// (aux weights, kept columns) onto the expanded/filtered edge set.
Rcpp::List build_geom_graph(const std::vector<double> &x1,
                             const std::vector<double> &y1,
                             const std::vector<double> &x2,
                             const std::vector<double> &y2,
                             const std::vector<double> &cost,
                             const std::vector<int> &direction,
                             int snap_digits);

#endif
