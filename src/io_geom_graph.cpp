#include "io_geom_graph.h"
#include <cmath>
#include <unordered_map>
#include <utility>
#include <fstream>
#include <iomanip>

namespace {

struct PairKey {
  long long a;
  long long b;
  bool operator==(const PairKey &other) const {
    return a == other.a && b == other.b;
  }
};

struct PairKeyHash {
  size_t operator()(const PairKey &k) const {
    size_t h = std::hash<long long>()(k.a);
    h ^= std::hash<long long>()(k.b) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    return h;
  }
};

// Core dedup/expansion algorithm, parameterized over how a node id is
// obtained/recorded (get_node_id) and how a finished edge is recorded
// (on_edge) -- shared by the in-memory List builder and the CSV streaming
// writer, so the endpoint-snapping/direction-expansion/self-loop-filtering
// logic exists in exactly one place either way. Unlike Module B, there is no
// "origin" index here: the CSV variant never carries aux (same restriction
// as the OSM CSV path), so on_edge doesn't need one.
template <typename GetNodeId, typename OnEdge>
void dedupe_and_split(const std::vector<double> &x1, const std::vector<double> &y1,
                       const std::vector<double> &x2, const std::vector<double> &y2,
                       const std::vector<double> &cost, const std::vector<int> &direction,
                       GetNodeId get_node_id, OnEdge on_edge){
  size_t n = x1.size();
  for (size_t i = 0; i < n; i++){
    int id1 = get_node_id(x1[i], y1[i]);
    int id2 = get_node_id(x2[i], y2[i]);
    if (id1 == id2) continue;

    if (direction[i] != 2) on_edge(id1, id2, cost[i], static_cast<int>(i));
    if (direction[i] != 1) on_edge(id2, id1, cost[i], static_cast<int>(i));
  }
}

}

Rcpp::List build_geom_graph(const std::vector<double> &x1,
                             const std::vector<double> &y1,
                             const std::vector<double> &x2,
                             const std::vector<double> &y2,
                             const std::vector<double> &cost,
                             const std::vector<int> &direction,
                             int snap_digits){

  double scale = std::pow(10.0, snap_digits);
  std::unordered_map<PairKey, int, PairKeyHash> node_map;
  node_map.reserve(x1.size() * 2);

  std::vector<double> node_x, node_y;
  auto get_node_id = [&](double x, double y) -> int {
    PairKey key{ static_cast<long long>(std::llround(x * scale)),
                 static_cast<long long>(std::llround(y * scale)) };
    auto it = node_map.find(key);
    if (it != node_map.end()) return it->second;
    int id = static_cast<int>(node_x.size());
    node_map.emplace(key, id);
    node_x.push_back(x);
    node_y.push_back(y);
    return id;
  };

  std::vector<int> from, to;
  std::vector<double> dist;
  std::vector<int> origin;
  size_t n = x1.size();
  from.reserve(n * 2); to.reserve(n * 2); dist.reserve(n * 2); origin.reserve(n * 2);

  auto on_edge = [&](int id1, int id2, double c, int orig){
    from.push_back(id1); to.push_back(id2); dist.push_back(c); origin.push_back(orig);
  };

  dedupe_and_split(x1, y1, x2, y2, cost, direction, get_node_id, on_edge);

  int nbnode = static_cast<int>(node_x.size());
  Rcpp::IntegerVector node_id(nbnode);
  for (int i = 0; i < nbnode; i++) node_id[i] = i;

  return Rcpp::List::create(
    Rcpp::Named("from") = Rcpp::wrap(from),
    Rcpp::Named("to") = Rcpp::wrap(to),
    Rcpp::Named("dist") = Rcpp::wrap(dist),
    Rcpp::Named("origin") = Rcpp::wrap(origin),
    Rcpp::Named("node_id") = node_id,
    Rcpp::Named("X") = Rcpp::wrap(node_x),
    Rcpp::Named("Y") = Rcpp::wrap(node_y)
  );
}

// [[Rcpp::export]]
Rcpp::List cpp_build_geom_graph(std::vector<double> x1, std::vector<double> y1,
                                 std::vector<double> x2, std::vector<double> y2,
                                 std::vector<double> cost, std::vector<int> direction,
                                 int snap_digits){
  return build_geom_graph(x1, y1, x2, y2, cost, direction, snap_digits);
}

// Streams "<edges_path>" (from,to,dist) and "<coords_path>" (node_id,X,Y)
// directly as each edge/node is produced, without ever holding the full
// edge/coordinate list in memory. Mirrors cpp_write_osm_graph_csv (Module
// B); the benefit is smaller here since Module A sources are typically far
// from continent-scale, but it's an option the user can still choose.
// [[Rcpp::export]]
Rcpp::List cpp_write_geom_graph_csv(std::vector<double> x1, std::vector<double> y1,
                                     std::vector<double> x2, std::vector<double> y2,
                                     std::vector<double> cost, std::vector<int> direction,
                                     int snap_digits,
                                     std::string edges_path, std::string coords_path){

  std::ofstream edges_f(edges_path.c_str());
  if (!edges_f) Rcpp::stop("Cannot open for writing: " + edges_path);
  edges_f << "from,to,dist\n" << std::setprecision(10);

  std::ofstream coords_f(coords_path.c_str());
  if (!coords_f) Rcpp::stop("Cannot open for writing: " + coords_path);
  coords_f << "node_id,X,Y\n" << std::setprecision(12);

  double scale = std::pow(10.0, snap_digits);
  std::unordered_map<PairKey, int, PairKeyHash> node_map;
  node_map.reserve(x1.size() * 2);

  auto get_node_id = [&](double x, double y) -> int {
    PairKey key{ static_cast<long long>(std::llround(x * scale)),
                 static_cast<long long>(std::llround(y * scale)) };
    auto it = node_map.find(key);
    if (it != node_map.end()) return it->second;
    int id = static_cast<int>(node_map.size());
    node_map.emplace(key, id);
    coords_f << id << "," << x << "," << y << "\n";
    return id;
  };

  auto on_edge = [&](int id1, int id2, double c, int /*orig*/){
    edges_f << id1 << "," << id2 << "," << c << "\n";
  };

  dedupe_and_split(x1, y1, x2, y2, cost, direction, get_node_id, on_edge);

  return Rcpp::List::create(Rcpp::Named("nbnode") = static_cast<int>(node_map.size()));
}
