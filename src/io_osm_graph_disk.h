#ifndef IO_OSM_GRAPH_DISK_H
#define IO_OSM_GRAPH_DISK_H

#include <Rcpp.h>
#include <string>
#include <vector>

// Disk-backed variant of build_osm_graph (io_osm_graph.h/.cpp): same
// tag-resolution and junction-splitting semantics, but bounded to roughly
// max_memory_mb of RAM regardless of the input file's total size, by
// spilling to temporary files (external sort + sequential merge-joins)
// instead of holding the referenced-node table / kept ways in RAM. Meant
// for .pbf files too large for the in-RAM path (io_osm_graph.h) on a given
// machine -- for anything that already fits comfortably (see the plan:
// country-scale extracts fit in a few GB), the in-RAM path is simpler and
// faster and remains the default.
//
// Always streams results to CSV (edges_path/coords_path): materializing the
// full result back into R would reintroduce the memory problem this exists
// to avoid, so unlike the in-RAM path there is no "graph"/"data.frame"
// output mode here -- the caller re-reads the CSV afterward if they know it
// will fit.
Rcpp::List build_osm_graph_disk(std::string path,
                                 std::vector<std::string> profile_highway,
                                 std::vector<int> profile_clazz,
                                 std::vector<double> profile_speed,
                                 bool use_maxspeed,
                                 bool assume_oneway_motorway,
                                 double max_memory_mb,
                                 std::string temp_dir,
                                 std::string edges_path,
                                 std::string coords_path,
                                 bool progress);

#endif
