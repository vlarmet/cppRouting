#include "io_geom_graph.h"
#include <cmath>
#include <unordered_map>
#include <utility>

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

  std::vector<double> node_x;
  std::vector<double> node_y;

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

  std::vector<int> from;
  std::vector<int> to;
  std::vector<double> dist;
  std::vector<int> origin;
  size_t n = x1.size();
  from.reserve(n * 2);
  to.reserve(n * 2);
  dist.reserve(n * 2);
  origin.reserve(n * 2);

  for (size_t i = 0; i < n; i++){
    int id1 = get_node_id(x1[i], y1[i]);
    int id2 = get_node_id(x2[i], y2[i]);
    if (id1 == id2) continue;

    if (direction[i] != 2){
      from.push_back(id1);
      to.push_back(id2);
      dist.push_back(cost[i]);
      origin.push_back(static_cast<int>(i));
    }
    if (direction[i] != 1){
      from.push_back(id2);
      to.push_back(id1);
      dist.push_back(cost[i]);
      origin.push_back(static_cast<int>(i));
    }
  }

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
