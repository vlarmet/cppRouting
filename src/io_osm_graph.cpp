#include "io_osm_graph.h"
#include "io_osmpbf.h"
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <cstdlib>
#include <cctype>

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

}

Rcpp::List build_osm_graph(std::string path,
                            std::vector<std::string> profile_highway,
                            std::vector<int> profile_clazz,
                            std::vector<double> profile_speed,
                            bool use_maxspeed,
                            bool assume_oneway_motorway,
                            bool progress){

  std::unordered_map<std::string, std::pair<int, double> > profile;
  for (size_t i = 0; i < profile_highway.size(); i++)
    profile[profile_highway[i]] = std::make_pair(profile_clazz[i], profile_speed[i]);

  // ---- Pass 1: filter ways, resolve tags, count node references ----
  if (progress) Rcpp::Rcout << "[cpp_extract] pass 1/3: reading ways..." << std::endl;

  std::vector<KeptWay> kept_ways;
  std::unordered_map<int64_t, uint32_t> ref_count;

  stream_osmpbf(path, std::function<void(const OsmNode&)>(),
    [&](const OsmWay &way){
      if (way.refs.size() < 2) return;
      const std::string *highway = find_tag(way.tags, "highway");
      if (!highway) return;
      auto it = profile.find(*highway);
      if (it == profile.end()) return;

      int clazz = it->second.first;
      double speed = it->second.second;

      bool allowed = true;
      const std::string *access = find_tag(way.tags, "access");
      if (access && (*access == "no" || *access == "private")) allowed = false;
      const std::string *motor_vehicle = find_tag(way.tags, "motor_vehicle");
      const std::string *motorcar = find_tag(way.tags, "motorcar");
      if (motor_vehicle && is_yesish(*motor_vehicle)) allowed = true;
      if (motorcar && is_yesish(*motorcar)) allowed = true;
      if (motor_vehicle && *motor_vehicle == "no") allowed = false;
      if (motorcar && *motorcar == "no") allowed = false;
      const std::string *vehicle = find_tag(way.tags, "vehicle");
      if (vehicle && *vehicle == "no") allowed = false;
      if (!allowed) return;

      if (use_maxspeed){
        const std::string *maxspeed = find_tag(way.tags, "maxspeed");
        double kmh;
        if (maxspeed && parse_maxspeed_kmh(*maxspeed, kmh) && kmh > 0) speed = kmh;
      }

      int direction = 0;
      const std::string *oneway = find_tag(way.tags, "oneway");
      if (oneway && (*oneway == "yes" || *oneway == "true" || *oneway == "1")) direction = 1;
      else if (oneway && (*oneway == "-1" || *oneway == "reverse")) direction = 2;
      else if (oneway && (*oneway == "no" || *oneway == "false" || *oneway == "0")) direction = 0;
      else if (assume_oneway_motorway && (*highway == "motorway" || *highway == "motorway_link")) direction = 1;

      kept_ways.push_back(KeptWay{way.refs, speed, direction, clazz});
      for (int64_t ref : way.refs) ref_count[ref]++;
    });

  // ---- Pass 2: collect coordinates for referenced nodes only ----
  if (progress) Rcpp::Rcout << "[cpp_extract] pass 2/3: reading node coordinates (" << ref_count.size() << " referenced nodes)..." << std::endl;

  std::unordered_map<int64_t, std::pair<float, float> > node_coord;
  node_coord.reserve(ref_count.size());

  stream_osmpbf(path,
    [&](const OsmNode &node){
      if (ref_count.find(node.id) != ref_count.end())
        node_coord[node.id] = std::make_pair(static_cast<float>(node.lat), static_cast<float>(node.lon));
    },
    std::function<void(const OsmWay&)>());

  // ---- Pass 3 (in-memory): split ways at junctions, build edges ----
  if (progress) Rcpp::Rcout << "[cpp_extract] pass 3/3: splitting ways at junctions and building edges..." << std::endl;

  std::unordered_map<int64_t, int> node_id_map;
  std::vector<double> node_lat_out, node_lon_out;
  auto get_node_id = [&](int64_t osm_id) -> int {
    auto it = node_id_map.find(osm_id);
    if (it != node_id_map.end()) return it->second;
    int id = static_cast<int>(node_lat_out.size());
    node_id_map.emplace(osm_id, id);
    auto coord = node_coord.at(osm_id);
    node_lat_out.push_back(coord.first);
    node_lon_out.push_back(coord.second);
    return id;
  };

  std::vector<int> from, to;
  std::vector<double> dist;
  std::vector<int> clazz_out;

  int n_boundary_ways_skipped = 0;
  for (const auto &way : kept_ways){
    // A regional extract can contain a way that crosses the extract's
    // boundary, referencing a node that was never included in the file at
    // all (not merely filtered out). Rather than crash on the missing
    // coordinate, drop the whole way -- it's a small, well-understood loss
    // limited to the extract's edge.
    bool all_present = true;
    for (int64_t ref : way.refs) if (node_coord.find(ref) == node_coord.end()) { all_present = false; break; }
    if (!all_present){ n_boundary_ways_skipped++; continue; }

    size_t n = way.refs.size();
    size_t seg_start = 0;
    for (size_t i = 1; i < n; i++){
      bool is_last = (i == n - 1);
      bool is_junction = ref_count.at(way.refs[i]) >= 2;
      if (!is_last && !is_junction) continue;

      int64_t from_osm = way.refs[seg_start];
      int64_t to_osm = way.refs[i];
      double length_m = 0;
      for (size_t j = seg_start; j < i; j++){
        auto c1 = node_coord.at(way.refs[j]);
        auto c2 = node_coord.at(way.refs[j + 1]);
        length_m += haversine_m(c1.first, c1.second, c2.first, c2.second);
      }

      int id1 = get_node_id(from_osm);
      int id2 = get_node_id(to_osm);
      if (id1 != id2){
        double cost = 60.0 * (length_m / 1000.0) / way.speed;
        if (way.direction != 2){ from.push_back(id1); to.push_back(id2); dist.push_back(cost); clazz_out.push_back(way.clazz); }
        if (way.direction != 1){ from.push_back(id2); to.push_back(id1); dist.push_back(cost); clazz_out.push_back(way.clazz); }
      }

      seg_start = i;
    }
  }

  int nbnode = static_cast<int>(node_lat_out.size());
  Rcpp::IntegerVector node_id(nbnode);
  for (int i = 0; i < nbnode; i++) node_id[i] = i;

  if (progress) Rcpp::Rcout << "[cpp_extract] done: " << nbnode << " nodes, " << from.size() << " edges ("
                            << n_boundary_ways_skipped << " way(s) skipped, crossing the extract's boundary)" << std::endl;

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
