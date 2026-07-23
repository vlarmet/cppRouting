#ifndef IO_OSM_TAGS_H
#define IO_OSM_TAGS_H

// Tag resolution (highway/access/maxspeed/oneway -> clazz/speed/direction)
// and the haversine distance helper, shared by the in-RAM OSM pipeline
// (io_osm_graph.cpp), the osmdata adapter (also io_osm_graph.cpp) and the
// disk-backed pipeline (io_osm_graph_disk.cpp) -- exactly one place for
// these rules so the three never drift apart.

#include "io_osmpbf.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <cstdlib>
#include <cctype>
#include <cmath>

namespace osm_tags {

inline const std::string* find_tag(const std::vector<OsmTag> &tags, const std::string &key){
  for (const auto &t : tags) if (t.key == key) return &t.value;
  return nullptr;
}

inline bool is_yesish(const std::string &v){
  return v == "yes" || v == "true" || v == "1" || v == "destination" || v == "permit";
}

// Parses a maxspeed value into km/h. Handles plain numbers (assumed km/h)
// and "<n> mph". Anything else (implicit zone speeds like "FR:urban",
// "none", "walk", ...) is left unresolved (caller keeps the class default).
inline bool parse_maxspeed_kmh(const std::string &v, double &out_kmh){
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

using Profile = std::unordered_map<std::string, std::pair<int, double> >;

inline Profile make_profile(const std::vector<std::string> &highway,
                             const std::vector<int> &clazz,
                             const std::vector<double> &speed){
  Profile profile;
  for (size_t i = 0; i < highway.size(); i++) profile[highway[i]] = std::make_pair(clazz[i], speed[i]);
  return profile;
}

// Resolves access tags for a given mode, following OSM's documented access
// hierarchy (access -> vehicle -> motor_vehicle -> mode-specific leaf tag;
// foot sits outside the vehicle branch entirely), least to most specific --
// a later, more specific check overrides an earlier one. The one exception
// is "car", whose check order (motor_vehicle/motorcar THEN vehicle) predates
// this generalization and is preserved exactly as-is for backward
// compatibility, even though it technically checks vehicle (less specific)
// after motorcar (more specific).
inline bool resolve_access(const std::vector<OsmTag> &tags, const std::string &mode){
  bool allowed = true;
  const std::string *access = find_tag(tags, "access");
  if (access && (*access == "no" || *access == "private")) allowed = false;

  if (mode == "foot"){
    const std::string *foot = find_tag(tags, "foot");
    if (foot && is_yesish(*foot)) allowed = true;
    if (foot && *foot == "no") allowed = false;
    return allowed;
  }

  if (mode == "bicycle"){
    const std::string *vehicle = find_tag(tags, "vehicle");
    if (vehicle && *vehicle == "no") allowed = false;
    const std::string *bicycle = find_tag(tags, "bicycle");
    if (bicycle && is_yesish(*bicycle)) allowed = true;
    if (bicycle && *bicycle == "no") allowed = false;
    return allowed;
  }

  // car and hgv both sit under motor_vehicle -- exact pre-existing car
  // order preserved (see comment above).
  const std::string *motor_vehicle = find_tag(tags, "motor_vehicle");
  const std::string *motorcar = find_tag(tags, "motorcar");
  if (motor_vehicle && is_yesish(*motor_vehicle)) allowed = true;
  if (motorcar && is_yesish(*motorcar)) allowed = true;
  if (motor_vehicle && *motor_vehicle == "no") allowed = false;
  if (motorcar && *motorcar == "no") allowed = false;
  const std::string *vehicle = find_tag(tags, "vehicle");
  if (vehicle && *vehicle == "no") allowed = false;
  if (mode != "hgv") return allowed;

  // hgv-specific tag checked last (most specific wins), e.g. hgv=destination
  // (is_yesish) re-authorizes, hgv=no overrides even motorcar=yes.
  const std::string *hgv = find_tag(tags, "hgv");
  if (hgv && is_yesish(*hgv)) allowed = true;
  if (hgv && *hgv == "no") allowed = false;
  return allowed;
}

// Resolves highway/access/maxspeed/oneway tags against the profile, for a
// given mode ("car" (default semantics), "foot", "bicycle", "hgv"). Returns
// false if the way should be dropped entirely.
inline bool resolve_way(const std::vector<OsmTag> &tags, const Profile &profile,
                         bool use_maxspeed, bool assume_oneway_motorway, const std::string &mode,
                         double &out_speed, int &out_direction, int &out_clazz){
  const std::string *highway = find_tag(tags, "highway");
  if (!highway) return false;
  auto it = profile.find(*highway);
  if (it == profile.end()) return false;

  out_clazz = it->second.first;
  out_speed = it->second.second;

  if (!resolve_access(tags, mode)) return false;

  // maxspeed is a legal speed limit for motor vehicles -- meaningless for a
  // pedestrian/cyclist's assumed pace, so only applied for car/hgv.
  if (use_maxspeed && (mode == "car" || mode == "hgv")){
    const std::string *maxspeed = find_tag(tags, "maxspeed");
    double kmh;
    if (maxspeed && parse_maxspeed_kmh(*maxspeed, kmh) && kmh > 0) out_speed = kmh;
  }

  // Pedestrians aren't meaningfully one-way: a rare explicit oneway:foot=yes
  // is not honored in this version (documented simplification).
  if (mode == "foot"){
    out_direction = 0;
    return true;
  }

  out_direction = 0;
  const std::string *oneway = find_tag(tags, "oneway");
  if (oneway && (*oneway == "yes" || *oneway == "true" || *oneway == "1")) out_direction = 1;
  else if (oneway && (*oneway == "-1" || *oneway == "reverse")) out_direction = 2;
  else if (oneway && (*oneway == "no" || *oneway == "false" || *oneway == "0")) out_direction = 0;
  else if (assume_oneway_motorway && (*highway == "motorway" || *highway == "motorway_link")) out_direction = 1;

  // Contraflow cycling: a one-way street can still be legally ridden the
  // "wrong" way by bicycle, tagged oneway:bicycle=no -- overrides whatever
  // the plain oneway tag said.
  if (mode == "bicycle"){
    const std::string *oneway_bicycle = find_tag(tags, "oneway:bicycle");
    if (oneway_bicycle && (*oneway_bicycle == "no" || *oneway_bicycle == "false" || *oneway_bicycle == "0")) out_direction = 0;
  }

  return true;
}

inline double haversine_m(double lat1, double lon1, double lat2, double lon2){
  const double R = 6371000.0;
  const double deg2rad = M_PI / 180.0;
  double dlat = (lat2 - lat1) * deg2rad;
  double dlon = (lon2 - lon1) * deg2rad;
  double a = std::sin(dlat / 2) * std::sin(dlat / 2) +
             std::cos(lat1 * deg2rad) * std::cos(lat2 * deg2rad) * std::sin(dlon / 2) * std::sin(dlon / 2);
  return 2 * R * std::atan2(std::sqrt(a), std::sqrt(1 - a));
}

}

#endif
