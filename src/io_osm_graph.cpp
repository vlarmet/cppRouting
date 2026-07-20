#include "io_osm_graph.h"
#include "io_osmpbf.h"
#include <cmath>
#include <limits>
#include <algorithm>
#include <unordered_map>
#include <cstdlib>
#include <cctype>
#include <fstream>
#include <iomanip>

namespace {

double haversine_m(double lat1, double lon1, double lat2, double lon2){
  const double R = 6371000.0;
  const double deg2rad = M_PI / 180.0;
  double dlat = (lat2 - lat1) * deg2rad;
  double dlon = (lon2 - lon1) * deg2rad;
  double a = std::sin(dlat / 2) * std::sin(dlat / 2) +
             std::cos(lat1 * deg2rad) * std::cos(lat2 * deg2rad) * std::sin(dlon / 2) * std::sin(dlon / 2);
  return 2 * R * std::atan2(std::sqrt(a), std::sqrt(1 - a));
}

const std::string* find_tag(const std::vector<OsmTag> &tags, const std::string &key){
  for (const auto &t : tags) if (t.key == key) return &t.value;
  return nullptr;
}

bool is_yesish(const std::string &v){
  return v == "yes" || v == "true" || v == "1" || v == "destination" || v == "permit";
}

// Parses a maxspeed value into km/h. Handles plain numbers (assumed km/h)
// and "<n> mph". Anything else (implicit zone speeds like "FR:urban",
// "none", "walk", ...) is left unresolved (caller keeps the class default).
bool parse_maxspeed_kmh(const std::string &v, double &out_kmh){
  size_t i = 0;
  while (i < v.size() && std::isspace(static_cast<unsigned char>(v[i]))) i++;
  size_t start = i;
  while (i < v.size() && (std::isdigit(static_cast<unsigned char>(v[i])) || v[i] == '.')) i++;
  if (i == start) return false;
  double num = std::atof(v.substr(start, i - start).c_str());
  while (i < v.size() && std::isspace(static_cast<unsigned char>(v[i]))) i++;
  std::string unit = v.substr(i);
  if (unit.empty()) { out_kmh = num; return true; }
  if (unit == "mph") { out_kmh = num * 1.609344; return true; }
  return false;
}

struct KeptWay {
  std::vector<int64_t> refs;
  double speed;
  int direction; // 0 = twoway, 1 = forward, 2 = backward
  int clazz;
};

// Resolves highway/access/maxspeed/oneway tags against the profile. Returns
// false if the way should be dropped entirely. Shared by both the .pbf path
// (tags come from the string table) and the osmdata adapter (tags come from
// sf attribute columns already as plain strings) -- same rules either way.
bool resolve_way(const std::vector<OsmTag> &tags,
                  const std::unordered_map<std::string, std::pair<int, double> > &profile,
                  bool use_maxspeed, bool assume_oneway_motorway,
                  double &out_speed, int &out_direction, int &out_clazz){
  const std::string *highway = find_tag(tags, "highway");
  if (!highway) return false;
  auto it = profile.find(*highway);
  if (it == profile.end()) return false;

  out_clazz = it->second.first;
  out_speed = it->second.second;

  bool allowed = true;
  const std::string *access = find_tag(tags, "access");
  if (access && (*access == "no" || *access == "private")) allowed = false;
  const std::string *motor_vehicle = find_tag(tags, "motor_vehicle");
  const std::string *motorcar = find_tag(tags, "motorcar");
  if (motor_vehicle && is_yesish(*motor_vehicle)) allowed = true;
  if (motorcar && is_yesish(*motorcar)) allowed = true;
  if (motor_vehicle && *motor_vehicle == "no") allowed = false;
  if (motorcar && *motorcar == "no") allowed = false;
  const std::string *vehicle = find_tag(tags, "vehicle");
  if (vehicle && *vehicle == "no") allowed = false;
  if (!allowed) return false;

  if (use_maxspeed){
    const std::string *maxspeed = find_tag(tags, "maxspeed");
    double kmh;
    if (maxspeed && parse_maxspeed_kmh(*maxspeed, kmh) && kmh > 0) out_speed = kmh;
  }

  out_direction = 0;
  const std::string *oneway = find_tag(tags, "oneway");
  if (oneway && (*oneway == "yes" || *oneway == "true" || *oneway == "1")) out_direction = 1;
  else if (oneway && (*oneway == "-1" || *oneway == "reverse")) out_direction = 2;
  else if (oneway && (*oneway == "no" || *oneway == "false" || *oneway == "0")) out_direction = 0;
  else if (assume_oneway_motorway && (*highway == "motorway" || *highway == "motorway_link")) out_direction = 1;

  return true;
}

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
      double length_m = 0;
      for (size_t j = seg_start; j < i; j++){
        auto c1 = nodes.coord(way.refs[j]);
        auto c2 = nodes.coord(way.refs[j + 1]);
        length_m += haversine_m(c1.first, c1.second, c2.first, c2.second);
      }

      int id1 = get_node_id(from_osm);
      int id2 = get_node_id(to_osm);
      if (id1 != id2){
        double cost = 60.0 * (length_m / 1000.0) / way.speed;
        if (way.direction != 2){ on_edge(id1, id2, cost, way.clazz); n_edges++; }
        if (way.direction != 1){ on_edge(id2, id1, cost, way.clazz); n_edges++; }
      }

      seg_start = i;
    }
  }

  return std::make_pair(n_edges, n_boundary_ways_skipped);
}

// Shared pass 3 (in-memory variant): returns from/to/dist/clazz plus a
// compact node_id/lon/lat table as an Rcpp::List. Used by the .pbf path and
// the osmdata adapter.
Rcpp::List build_edges_from_ways(const std::vector<KeptWay> &kept_ways,
                                  const NodeTable &nodes,
                                  bool progress){
  std::unordered_map<int64_t, int> node_id_map;
  std::vector<double> node_lat_out, node_lon_out;
  auto get_node_id = [&](int64_t osm_id) -> int {
    auto it = node_id_map.find(osm_id);
    if (it != node_id_map.end()) return it->second;
    int id = static_cast<int>(node_lat_out.size());
    node_id_map.emplace(osm_id, id);
    auto coord = nodes.coord(osm_id);
    node_lat_out.push_back(coord.first);
    node_lon_out.push_back(coord.second);
    return id;
  };

  std::vector<int> from, to;
  std::vector<double> dist;
  std::vector<int> clazz_out;
  auto on_edge = [&](int id1, int id2, double cost, int clazz){
    from.push_back(id1); to.push_back(id2); dist.push_back(cost); clazz_out.push_back(clazz);
  };

  auto counts = split_ways_into_edges(kept_ways, nodes, get_node_id, on_edge);

  int nbnode = static_cast<int>(node_lat_out.size());
  Rcpp::IntegerVector node_id(nbnode);
  for (int i = 0; i < nbnode; i++) node_id[i] = i;

  if (progress) Rcpp::Rcout << "[cpp_extract] done: " << nbnode << " nodes, " << counts.first << " edges ("
                            << counts.second << " way(s) skipped, crossing the extract's boundary)" << std::endl;

  return Rcpp::List::create(
    Rcpp::Named("from") = Rcpp::wrap(from),
    Rcpp::Named("to") = Rcpp::wrap(to),
    Rcpp::Named("dist") = Rcpp::wrap(dist),
    Rcpp::Named("clazz") = Rcpp::wrap(clazz_out),
    Rcpp::Named("node_id") = node_id,
    Rcpp::Named("lon") = Rcpp::wrap(node_lon_out),
    Rcpp::Named("lat") = Rcpp::wrap(node_lat_out)
  );
}

// Shared pass 3 (CSV streaming variant): writes "<edges_path>" (from,to,dist)
// and "<coords_path>" (node_id,X,Y) directly as each edge/node is produced,
// without ever holding the full edge/coordinate list in memory -- the point
// being to keep peak RAM bounded by pass 1/2's structures alone, for
// continent-scale files processed with output="csv".
Rcpp::List write_edges_from_ways_csv(const std::vector<KeptWay> &kept_ways,
                                      const NodeTable &nodes,
                                      const std::string &edges_path,
                                      const std::string &coords_path,
                                      bool progress){
  std::ofstream edges_f(edges_path.c_str());
  if (!edges_f) Rcpp::stop("Cannot open for writing: " + edges_path);
  edges_f << "from,to,dist\n" << std::setprecision(10);

  std::ofstream coords_f(coords_path.c_str());
  if (!coords_f) Rcpp::stop("Cannot open for writing: " + coords_path);
  coords_f << "node_id,X,Y\n" << std::setprecision(9);

  std::unordered_map<int64_t, int> node_id_map;
  auto get_node_id = [&](int64_t osm_id) -> int {
    auto it = node_id_map.find(osm_id);
    if (it != node_id_map.end()) return it->second;
    int id = static_cast<int>(node_id_map.size());
    node_id_map.emplace(osm_id, id);
    auto coord = nodes.coord(osm_id);
    coords_f << id << "," << coord.second << "," << coord.first << "\n"; // X=lon, Y=lat
    return id;
  };

  auto on_edge = [&](int id1, int id2, double cost, int /*clazz*/){
    edges_f << id1 << "," << id2 << "," << cost << "\n";
  };

  auto counts = split_ways_into_edges(kept_ways, nodes, get_node_id, on_edge);
  int nbnode = static_cast<int>(node_id_map.size());

  if (progress) Rcpp::Rcout << "[cpp_extract] done: " << nbnode << " nodes, " << counts.first << " edges ("
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
  std::unordered_map<std::string, std::pair<int, double> > profile;
  for (size_t i = 0; i < profile_highway.size(); i++)
    profile[profile_highway[i]] = std::make_pair(profile_clazz[i], profile_speed[i]);

  if (progress) Rcpp::Rcout << "[cpp_extract] pass 1/3: reading ways..." << std::endl;

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

  if (progress) Rcpp::Rcout << "[cpp_extract] pass 2/3: reading node coordinates (" << nodes.size() << " referenced nodes)..." << std::endl;

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

  if (progress) Rcpp::Rcout << "[cpp_extract] pass 3/3: splitting ways at junctions and building edges..." << std::endl;
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
                                     bool progress){
  std::vector<KeptWay> kept_ways;
  NodeTable nodes;
  run_passes_1_2(path, profile_highway, profile_clazz, profile_speed,
                 use_maxspeed, assume_oneway_motorway, progress, kept_ways, nodes);

  if (progress) Rcpp::Rcout << "[cpp_extract] pass 3/3: splitting ways at junctions and streaming edges to CSV..." << std::endl;
  return write_edges_from_ways_csv(kept_ways, nodes, edges_path, coords_path, progress);
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

  std::unordered_map<std::string, std::pair<int, double> > profile;
  for (size_t i = 0; i < profile_highway.size(); i++)
    profile[profile_highway[i]] = std::make_pair(profile_clazz[i], profile_speed[i]);

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
