#ifndef IO_OSM_GRAPH_H
#define IO_OSM_GRAPH_H

#include <Rcpp.h>
#include <string>
#include <vector>

// Builds a from/to/dist/coords edge table from an .osm.pbf file (Module B):
// filters ways using the given highway->clazz/speed profile, resolves
// maxspeed/oneway/access tags, and detects true junctions (a node
// referenced by >= 2 kept ways, or a way endpoint) to split ways into
// edges -- unlike Module A, topology comes from OSM's own stable node ids,
// not coordinate matching, since a way can pass through a junction in its
// middle.
Rcpp::List build_osm_graph(std::string path,
                            std::vector<std::string> profile_highway,
                            std::vector<int> profile_clazz,
                            std::vector<double> profile_speed,
                            bool use_maxspeed,
                            bool assume_oneway_motorway,
                            bool progress);

#endif
