#include "io_osm_graph.h"
#include "io_osmpbf.h"
#include "io_osm_tags.h"
#include <cmath>
#include <limits>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <cstdlib>
#include <cctype>
#include <fstream>
#include <iomanip>

using osm_tags::haversine_m;
using osm_tags::resolve_way;
using osm_tags::make_profile;

namespace {

struct KeptWay {
  std::vector<int64_t> refs;
  double speed;
  int direction; // 0 = twoway, 1 = forward, 2 = backward
  int clazz;
};

// Sorted, deduplicated node table: id + coordinate + reference count, one
// contiguous array with O(log n) binary-search lookup. Replaces what used
// to be two separate unordered_maps (ref_count + node_coord). Those had
// per-entry hash-map overhead (chaining pointers, bucket bookkeeping) on top
// of the 16 bytes of actual payload -- for the tens of millions of nodes a
// country/continent extract references, that overhead was the dominant cost.
// A sorted vector needs exactly the payload, no more, and -- unlike a hash
// map -- is a single contiguous block that could later be backed by an
// mmap'd file for a true disk mode, if a file ever gets too big for RAM.
class NodeTable {
public:
  // Consumes raw_refs (sorts it, then frees it): one entry per (way, node)
  // occurrence, duplicates expected. Collapsed into one entry per unique id,
  // with ref_count = number of occurrences. Coordinates start unset (NaN)
  // until set_coord() is called for each one during pass 2.
  void build_from_refs(std::vector<int64_t> &raw_refs){
    std::sort(raw_refs.begin(), raw_refs.end());
    nodes_.clear();
    nodes_.reserve(raw_refs.size());
    size_t i = 0;
    while (i < raw_refs.size()){
      size_t j = i;
      while (j < raw_refs.size() && raw_refs[j] == raw_refs[i]) j++;
      nodes_.push_back(Entry{raw_refs[i], kUnset, kUnset, static_cast<uint32_t>(j - i)});
      i = j;
    }
    std::vector<int64_t>().swap(raw_refs);
  }

  void set_coord(int64_t id, float lat, float lon){
    Entry *e = find_mut(id);
    if (e) { e->lat = lat; e->lon = lon; }
  }

  bool has_coord(int64_t id) const {
    const Entry *e = find(id);
    return e != nullptr && e->lat != kUnset;
  }

  uint32_t ref_count(int64_t id) const {
    const Entry *e = find(id);
    return e ? e->ref_count : 0;
  }

  // Caller must have already checked has_coord(id).
  std::pair<float, float> coord(int64_t id) const {
    const Entry *e = find(id);
    return std::make_pair(e->lat, e->lon);
  }

  size_t size() const { return nodes_.size(); }

private:
  static constexpr float kUnset = std::numeric_limits<float>::lowest();

  struct Entry { int64_t id; float lat; float lon; uint32_t ref_count; };

  const Entry* find(int64_t id) const {
    auto it = std::lower_bound(nodes_.begin(), nodes_.end(), id,
                                [](const Entry &e, int64_t v){ return e.id < v; });
    if (it != nodes_.end() && it->id == id) return &(*it);
    return nullptr;
  }
  Entry* find_mut(int64_t id){
    auto it = std::lower_bound(nodes_.begin(), nodes_.end(), id,
                                [](const Entry &e, int64_t v){ return e.id < v; });
    if (it != nodes_.end() && it->id == id) return &(*it);
    return nullptr;
  }

  std::vector<Entry> nodes_;
};

// Core junction-splitting algorithm, parameterized over how a node id is
// obtained/recorded (get_node_id) and how a finished edge is recorded
// (on_edge) -- shared by the in-memory List builder (nodes/edges pushed
// into vectors) and the CSV streaming writer (nodes/edges written straight
// to a file the moment they're produced, so neither ever needs to be held
// as one big in-memory vector). The junction rule itself (a node referenced
// by >= 2 kept ways, or a way's own endpoint) and the haversine length
// computation exist in exactly this one place either way.
template <typename GetNodeId, typename OnEdge>
std::pair<long long, int> split_ways_into_edges(const std::vector<KeptWay> &kept_ways,
                                                 const NodeTable &nodes,
                                                 GetNodeId get_node_id, OnEdge on_edge){
  long long n_edges = 0;
  int n_boundary_ways_skipped = 0;

  for (const auto &way : kept_ways){
    // A regional extract can contain a way that crosses the extract's
    // boundary, referencing a node that was never included in the file at
    // all (not merely filtered out). Rather than crash on the missing
    // coordinate, drop the whole way -- it's a small, well-understood loss
    // limited to the extract's edge.
    bool all_present = true;
    for (int64_t ref : way.refs) if (!nodes.has_coord(ref)) { all_present = false; break; }
    if (!all_present){ n_boundary_ways_skipped++; continue; }

    size_t n = way.refs.size();
    size_t seg_start = 0;
    for (size_t i = 1; i < n; i++){
      bool is_last = (i == n - 1);
      bool is_junction = nodes.ref_count(way.refs[i]) >= 2;
      if (!is_last && !is_junction) continue;

      int64_t from_osm = way.refs[seg_start];
      int64_t to_osm = way.refs[i];

      // Checked on the OSM ids themselves (equivalent to id1 == id2 after
      // get_node_id, since that mapping is a bijection) so a self-loop never
      // registers its endpoint at all -- matches the disk pipeline's
      // behaviour (io_osm_graph_disk.cpp), which never sees these ids
      // either since it only ever adds a coord when from_osm != to_osm.
      if (from_osm != to_osm){
        double length_m = 0;
        for (size_t j = seg_start; j < i; j++){
          auto c1 = nodes.coord(way.refs[j]);
          auto c2 = nodes.coord(way.refs[j + 1]);
          length_m += haversine_m(c1.first, c1.second, c2.first, c2.second);
        }

        auto id1 = get_node_id(from_osm);
        auto id2 = get_node_id(to_osm);
        double cost = 60.0 * (length_m / 1000.0) / way.speed;
        double length_km = length_m / 1000.0;
        if (way.direction != 2){ on_edge(id1, id2, cost, way.clazz, length_km); n_edges++; }
        if (way.direction != 1){ on_edge(id2, id1, cost, way.clazz, length_km); n_edges++; }
      }

      seg_start = i;
    }
  }

  return std::make_pair(n_edges, n_boundary_ways_skipped);
}

// Shared pass 3 (in-memory variant): returns from/to/dist/clazz plus a
// node_id/lon/lat table as an Rcpp::List. Used by the .pbf path and the
// osmdata adapter. from/to/node_id carry the real OSM node id (as double --
// current OSM ids already exceed 2^31, so plain `int`/IntegerVector would
// silently overflow; double is exact up to 2^53, comfortably enough
// headroom) rather than a synthetic dense index: OSM ids are already stable,
// unique identifiers, so there's no reason to discard them the way Module A
// (shapefile/sf, which has no such native id) has to.
Rcpp::List build_edges_from_ways(const std::vector<KeptWay> &kept_ways,
                                  const NodeTable &nodes,
                                  bool progress){
  std::unordered_set<int64_t> seen;
  std::vector<double> node_id_out, node_lat_out, node_lon_out;
  auto get_node_id = [&](int64_t osm_id) -> int64_t {
    if (seen.insert(osm_id).second) {
      auto coord = nodes.coord(osm_id);
      node_id_out.push_back(static_cast<double>(osm_id));
      node_lat_out.push_back(coord.first);
      node_lon_out.push_back(coord.second);
    }
    return osm_id;
  };

  std::vector<double> from, to;
  std::vector<double> dist;
  std::vector<int> clazz_out;
  std::vector<double> length_out;
  auto on_edge = [&](int64_t id1, int64_t id2, double cost, int clazz, double length_km){
    from.push_back(static_cast<double>(id1)); to.push_back(static_cast<double>(id2));
    dist.push_back(cost); clazz_out.push_back(clazz);
    length_out.push_back(length_km);
  };

  auto counts = split_ways_into_edges(kept_ways, nodes, get_node_id, on_edge);

  int nbnode = static_cast<int>(node_id_out.size());

  if (progress) Rcpp::Rcout << "done: " << nbnode << " nodes, " << counts.first << " edges ("
                            << counts.second << " way(s) skipped, crossing the extract's boundary)" << std::endl;

  return Rcpp::List::create(
    Rcpp::Named("from") = Rcpp::wrap(from),
    Rcpp::Named("to") = Rcpp::wrap(to),
    Rcpp::Named("length") = Rcpp::wrap(length_out),
    Rcpp::Named("dist") = Rcpp::wrap(dist),
    Rcpp::Named("clazz") = Rcpp::wrap(clazz_out),
    Rcpp::Named("node_id") = Rcpp::wrap(node_id_out),
    Rcpp::Named("lon") = Rcpp::wrap(node_lon_out),
    Rcpp::Named("lat") = Rcpp::wrap(node_lat_out)
  );
}

// Shared pass 3 (CSV streaming variant): writes "<edges_path>" (from,to,dist),
// if want_coords, "<coords_path>" (node_id,X,Y), and if keep_clazz/keep_length,
// "<aux_path>" (clazz and/or length) directly as each edge/node is produced,
// without ever holding the full edge/coordinate/aux list in memory -- the
// point being to keep peak RAM bounded by pass 1/2's structures alone, for
// continent-scale files processed with output="csv". aux_path mirrors
// data.frame mode's separate `aux` table exactly: same columns, same values,
// just written straight to disk instead of returned as an R object. from/to/
// node_id carry the real OSM node id (matching the max_memory disk pipeline
// in io_osm_graph_disk.cpp, and build_edges_from_ways() above) rather than a
// synthetic dense index. When want_coords is false, coords_path is never
// created -- the `seen` set is still needed regardless (it's what dedups a
// node written once per way that uses it as a split-point endpoint), only
// the coordinate write itself is skipped. Likewise aux_path is never created
// unless at least one of keep_clazz/keep_length is requested.
Rcpp::List write_edges_from_ways_csv(const std::vector<KeptWay> &kept_ways,
                                      const NodeTable &nodes,
                                      const std::string &edges_path,
                                      const std::string &coords_path,
                                      const std::string &aux_path,
                                      bool want_coords,
                                      bool keep_clazz,
                                      bool keep_length,
                                      bool progress){
  std::ofstream edges_f(edges_path.c_str());
  if (!edges_f) Rcpp::stop("Cannot open for writing: " + edges_path);
  edges_f << "from,to,dist\n" << std::setprecision(10);

  std::ofstream coords_f;
  if (want_coords) {
    coords_f.open(coords_path.c_str());
    if (!coords_f) Rcpp::stop("Cannot open for writing: " + coords_path);
    coords_f << "node_id,X,Y\n" << std::setprecision(9);
  }

  bool want_aux = keep_clazz || keep_length;
  std::ofstream aux_f;
  if (want_aux) {
    aux_f.open(aux_path.c_str());
    if (!aux_f) Rcpp::stop("Cannot open for writing: " + aux_path);
    aux_f << std::setprecision(10);
    if (keep_clazz && keep_length) aux_f << "clazz,length\n";
    else if (keep_clazz) aux_f << "clazz\n";
    else aux_f << "length\n";
  }

  // Real OSM node ids already exceed 2^31: written as plain digits, they'd
  // look like a valid (if large) integer to data.table::fread(), which then
  // silently promotes the column to bit64::integer64 -- and without the
  // (Suggests-only, not required by this package) bit64 package installed,
  // that class holds outright garbage (an int64 bit pattern reinterpreted as
  // a plain double), not merely a display quirk. A trailing ".0" makes the
  // column unambiguously a double on read, which is exact for any id up to
  // 2^53 -- comfortably past any real OSM id.
  std::unordered_set<int64_t> seen;
  auto get_node_id = [&](int64_t osm_id) -> int64_t {
    if (seen.insert(osm_id).second && want_coords) {
      auto coord = nodes.coord(osm_id);
      coords_f << osm_id << ".0," << coord.second << "," << coord.first << "\n"; // X=lon, Y=lat
    }
    return osm_id;
  };

  auto on_edge = [&](int64_t id1, int64_t id2, double cost, int clazz, double length_km){
    edges_f << id1 << ".0," << id2 << ".0," << cost << "\n";
    if (want_aux) {
      if (keep_clazz && keep_length) aux_f << clazz << "," << length_km << "\n";
      else if (keep_clazz) aux_f << clazz << "\n";
      else aux_f << length_km << "\n";
    }
  };

  auto counts = split_ways_into_edges(kept_ways, nodes, get_node_id, on_edge);
  int nbnode = static_cast<int>(seen.size());

  if (progress) Rcpp::Rcout << "done: " << nbnode << " nodes, " << counts.first << " edges ("
                            << counts.second << " way(s) skipped, crossing the extract's boundary)" << std::endl;

  return Rcpp::List::create(Rcpp::Named("nbnode") = nbnode, Rcpp::Named("nbedge") = static_cast<double>(counts.first));
}

// Passes 1 and 2 (filter+resolve ways, then collect coordinates for the
// nodes they reference) are identical whether pass 3 builds an in-memory
// List or streams straight to CSV -- shared here so the two only ever
// diverge at pass 3. Pass 1 collects raw (duplicated) node references
// instead of counting them incrementally into a map, so the NodeTable can
// be built once as a sorted, deduplicated block (see build_from_refs).
void run_passes_1_2(const std::string &path,
                     const std::vector<std::string> &profile_highway,
                     const std::vector<int> &profile_clazz,
                     const std::vector<double> &profile_speed,
                     bool use_maxspeed, bool assume_oneway_motorway, bool progress,
                     std::vector<KeptWay> &kept_ways,
                     NodeTable &nodes){
  osm_tags::Profile profile = make_profile(profile_highway, profile_clazz, profile_speed);

  if (progress) Rcpp::Rcout << "pass 1/3: reading ways..." << std::endl;

  std::vector<int64_t> raw_refs;
  stream_osmpbf(path, std::function<void(const OsmNode&)>(),
    [&](const OsmWay &way){
      if (way.refs.size() < 2) return;
      double speed; int direction; int clazz;
      if (!resolve_way(way.tags, profile, use_maxspeed, assume_oneway_motorway, speed, direction, clazz)) return;
      kept_ways.push_back(KeptWay{way.refs, speed, direction, clazz});
      raw_refs.insert(raw_refs.end(), way.refs.begin(), way.refs.end());
    });

  nodes.build_from_refs(raw_refs);

  if (progress) Rcpp::Rcout << "pass 2/3: reading node coordinates (" << nodes.size() << " referenced nodes)..." << std::endl;

  stream_osmpbf(path,
    [&](const OsmNode &node){ nodes.set_coord(node.id, static_cast<float>(node.lat), static_cast<float>(node.lon)); },
    std::function<void(const OsmWay&)>());
}

}

Rcpp::List build_osm_graph(std::string path,
                            std::vector<std::string> profile_highway,
                            std::vector<int> profile_clazz,
                            std::vector<double> profile_speed,
                            bool use_maxspeed,
                            bool assume_oneway_motorway,
                            bool progress){
  std::vector<KeptWay> kept_ways;
  NodeTable nodes;
  run_passes_1_2(path, profile_highway, profile_clazz, profile_speed,
                 use_maxspeed, assume_oneway_motorway, progress, kept_ways, nodes);

  if (progress) Rcpp::Rcout << "pass 3/3: splitting ways at junctions and building edges..." << std::endl;
  return build_edges_from_ways(kept_ways, nodes, progress);
}

// [[Rcpp::export]]
Rcpp::List cpp_write_osm_graph_csv(std::string path,
                                     std::vector<std::string> profile_highway,
                                     std::vector<int> profile_clazz,
                                     std::vector<double> profile_speed,
                                     bool use_maxspeed,
                                     bool assume_oneway_motorway,
                                     std::string edges_path,
                                     std::string coords_path,
                                     std::string aux_path,
                                     bool want_coords,
                                     bool keep_clazz,
                                     bool keep_length,
                                     bool progress){
  std::vector<KeptWay> kept_ways;
  NodeTable nodes;
  run_passes_1_2(path, profile_highway, profile_clazz, profile_speed,
                 use_maxspeed, assume_oneway_motorway, progress, kept_ways, nodes);

  if (progress) Rcpp::Rcout << "pass 3/3: splitting ways at junctions and streaming edges to CSV..." << std::endl;
  return write_edges_from_ways_csv(kept_ways, nodes, edges_path, coords_path, aux_path,
                                    want_coords, keep_clazz, keep_length, progress);
}

// [[Rcpp::export]]
Rcpp::List cpp_build_osm_graph(std::string path,
                                 std::vector<std::string> profile_highway,
                                 std::vector<int> profile_clazz,
                                 std::vector<double> profile_speed,
                                 bool use_maxspeed,
                                 bool assume_oneway_motorway,
                                 bool progress){
  return build_osm_graph(path, profile_highway, profile_clazz, profile_speed,
                          use_maxspeed, assume_oneway_motorway, progress);
}

// osmdata adapter entry point: way_refs/way_tag_keys/way_tag_vals are
// parallel R lists (one element per way -- refs a numeric vector of OSM
// node ids reconstructed by the R side via coordinate-matching against
// osm_points, tag_keys/tag_vals character vectors of equal length). Node
// coordinates are provided directly (from osm_points), not streamed from a
// file, since osmdata results are already fully decoded in memory.
// [[Rcpp::export]]
Rcpp::List cpp_build_osm_graph_from_data(Rcpp::List way_refs,
                                          Rcpp::List way_tag_keys,
                                          Rcpp::List way_tag_vals,
                                          std::vector<double> node_ids,
                                          std::vector<double> node_lats,
                                          std::vector<double> node_lons,
                                          std::vector<std::string> profile_highway,
                                          std::vector<int> profile_clazz,
                                          std::vector<double> profile_speed,
                                          bool use_maxspeed,
                                          bool assume_oneway_motorway){

  osm_tags::Profile profile = make_profile(profile_highway, profile_clazz, profile_speed);

  std::vector<KeptWay> kept_ways;
  std::vector<int64_t> raw_refs;

  int n_ways = way_refs.size();
  for (int w = 0; w < n_ways; w++){
    std::vector<double> refs_d = Rcpp::as<std::vector<double> >(way_refs[w]);
    if (refs_d.size() < 2) continue;
    std::vector<std::string> keys = Rcpp::as<std::vector<std::string> >(way_tag_keys[w]);
    std::vector<std::string> vals = Rcpp::as<std::vector<std::string> >(way_tag_vals[w]);

    std::vector<OsmTag> tags;
    tags.reserve(keys.size());
    for (size_t k = 0; k < keys.size() && k < vals.size(); k++) tags.push_back(OsmTag{keys[k], vals[k]});

    double speed; int direction; int clazz;
    if (!resolve_way(tags, profile, use_maxspeed, assume_oneway_motorway, speed, direction, clazz)) continue;

    std::vector<int64_t> refs(refs_d.size());
    for (size_t i = 0; i < refs_d.size(); i++) refs[i] = static_cast<int64_t>(refs_d[i]);

    kept_ways.push_back(KeptWay{refs, speed, direction, clazz});
    raw_refs.insert(raw_refs.end(), refs.begin(), refs.end());
  }

  NodeTable nodes;
  nodes.build_from_refs(raw_refs);
  for (size_t i = 0; i < node_ids.size(); i++)
    nodes.set_coord(static_cast<int64_t>(node_ids[i]), static_cast<float>(node_lats[i]), static_cast<float>(node_lons[i]));

  return build_edges_from_ways(kept_ways, nodes, false);
}
