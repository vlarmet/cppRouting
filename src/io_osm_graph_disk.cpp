#include "io_osm_graph_disk.h"
#include "io_osmpbf.h"
#include "io_osm_tags.h"
#include "external_sort.h"
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cstdint>
#include <functional>
#include <memory>

namespace {

// ---- Fixed-size POD records for the external sorts / merge-joins ----

// (node_id, way index, position within the way) -- sorted by node_id to be
// merge-joined against the node coordinate table.
// explicit trailing pad: the record is serialized as raw bytes (ExternalSorter),
// so implicit alignment padding would carry uninitialised bytes to disk
struct IdPosRecord { int64_t node_id; int64_t way_idx; int32_t position; int32_t pad; };
struct CompareIdPosById {
  bool operator()(const IdPosRecord &a, const IdPosRecord &b) const { return a.node_id < b.node_id; }
};

// (way index, position, coordinate, ref_count) -- sorted by (way_idx,
// position) to reconstruct each way's node sequence in original order.
struct WayPosRecord { int64_t way_idx; int32_t position; float lat; float lon; uint32_t ref_count; };
struct CompareWayPos {
  bool operator()(const WayPosRecord &a, const WayPosRecord &b) const {
    if (a.way_idx != b.way_idx) return a.way_idx < b.way_idx;
    return a.position < b.position;
  }
};

// (node_id, coordinate) for the final graph-node table -- collected once
// per way-endpoint occurrence (so duplicates across ways are expected) and
// deduplicated in one pass at the very end.
struct CoordRecord { int64_t node_id; float lat; float lon; };
struct CompareCoordById {
  bool operator()(const CoordRecord &a, const CoordRecord &b) const { return a.node_id < b.node_id; }
};

// ---- Sequential (non-sorted) on-disk storage for kept ways ----
// Ways are read/written strictly in the order they were kept in pass 1;
// that position IS the way_idx used throughout the pipeline, so no index
// needs to be stored explicitly.

class WaysFileWriter {
public:
  explicit WaysFileWriter(const std::string &path) : f_(path.c_str(), std::ios::binary) {
    if (!f_) Rcpp::stop("Cannot create temp file: " + path);
  }
  void write_way(double speed, int direction, int clazz, const std::vector<int64_t> &refs){
    int32_t dir32 = direction, clazz32 = clazz, n = static_cast<int32_t>(refs.size());
    f_.write(reinterpret_cast<const char*>(&speed), sizeof(double));
    f_.write(reinterpret_cast<const char*>(&dir32), sizeof(int32_t));
    f_.write(reinterpret_cast<const char*>(&clazz32), sizeof(int32_t));
    f_.write(reinterpret_cast<const char*>(&n), sizeof(int32_t));
    f_.write(reinterpret_cast<const char*>(refs.data()), static_cast<std::streamsize>(refs.size() * sizeof(int64_t)));
  }
private:
  std::ofstream f_;
};

class WaysFileReader {
public:
  explicit WaysFileReader(const std::string &path) : f_(path.c_str(), std::ios::binary) {
    if (!f_) Rcpp::stop("Cannot open temp file: " + path);
  }
  bool read_way(double &speed, int &direction, int &clazz, std::vector<int64_t> &refs){
    if (!f_.read(reinterpret_cast<char*>(&speed), sizeof(double))) return false;
    int32_t dir32, clazz32, n;
    f_.read(reinterpret_cast<char*>(&dir32), sizeof(int32_t));
    f_.read(reinterpret_cast<char*>(&clazz32), sizeof(int32_t));
    f_.read(reinterpret_cast<char*>(&n), sizeof(int32_t));
    direction = dir32; clazz = clazz32;
    refs.resize(static_cast<size_t>(n));
    f_.read(reinterpret_cast<char*>(refs.data()), static_cast<std::streamsize>(refs.size() * sizeof(int64_t)));
    return true;
  }
private:
  std::ifstream f_;
};

}

Rcpp::List build_osm_graph_disk(std::string path,
                                 std::vector<std::string> profile_highway,
                                 std::vector<int> profile_clazz,
                                 std::vector<double> profile_speed,
                                 bool use_maxspeed,
                                 bool assume_oneway_motorway,
                                 std::string mode,
                                 double max_memory_mb,
                                 std::string temp_dir,
                                 std::string edges_path,
                                 std::string coords_path,
                                 std::string aux_path,
                                 bool want_coords,
                                 bool keep_clazz,
                                 bool keep_length,
                                 bool progress){

  osm_tags::Profile profile = osm_tags::make_profile(profile_highway, profile_clazz, profile_speed);

  size_t bytes_budget = static_cast<size_t>(max_memory_mb * 1024.0 * 1024.0);
  size_t ref_chunk    = std::max<size_t>(1024, bytes_budget / sizeof(int64_t));
  size_t idpos_chunk  = std::max<size_t>(1024, bytes_budget / sizeof(IdPosRecord));
  size_t waypos_chunk = std::max<size_t>(1024, bytes_budget / sizeof(WayPosRecord));
  size_t coord_chunk  = std::max<size_t>(1024, bytes_budget / sizeof(CoordRecord));

  std::string ways_path       = temp_dir + "/cppRouting_osm_ways.bin";
  std::string refcount_path   = temp_dir + "/cppRouting_osm_refcount.bin";
  std::string coords_bin_path = temp_dir + "/cppRouting_osm_coords.bin";

  // ---- Pass 1/5: filter ways, stream them to disk, collect raw node refs ----
  if (progress) Rcpp::Rcout << "pass 1/5: reading ways..." << std::endl;

  ExternalSorter<int64_t> ref_sorter(ref_chunk, temp_dir, "cppRouting_osm_refs");
  long long n_kept_ways = 0;
  {
    WaysFileWriter ways_writer(ways_path);
    stream_osmpbf(path, std::function<void(const OsmNode&)>(),
      [&](const OsmWay &way){
        if (way.refs.size() < 2) return;
        double speed; int direction; int clazz;
        if (!osm_tags::resolve_way(way.tags, profile, use_maxspeed, assume_oneway_motorway, mode, speed, direction, clazz)) return;
        ways_writer.write_way(speed, direction, clazz, way.refs);
        for (int64_t ref : way.refs) ref_sorter.add(ref);
        n_kept_ways++;
      });
  }
  ref_sorter.finish();

  // ---- Pass 2/5: collapse the sorted (duplicated) ref stream into (id, ref_count), sorted by id ----
  if (progress) Rcpp::Rcout << "pass 2/5: deduplicating referenced nodes..." << std::endl;

  long long n_distinct_refs = 0;
  {
    std::ofstream rc_f(refcount_path.c_str(), std::ios::binary);
    if (!rc_f) Rcpp::stop("Cannot create temp file: " + refcount_path);
    auto it = ref_sorter.iterator();
    if (it.has_next()){
      int64_t cur = it.next();
      uint32_t count = 1;
      auto flush_current = [&](){
        rc_f.write(reinterpret_cast<const char*>(&cur), sizeof(int64_t));
        rc_f.write(reinterpret_cast<const char*>(&count), sizeof(uint32_t));
        n_distinct_refs++;
      };
      while (it.has_next()){
        int64_t v = it.next();
        if (v == cur) { count++; }
        else { flush_current(); cur = v; count = 1; }
      }
      flush_current();
    }
  }

  // ---- Pass 3/5: merge-join the (already ~sorted) pbf node stream against node_refcount, -> node_coords.bin (sorted by id) ----
  if (progress) Rcpp::Rcout << "pass 3/5: reading node coordinates (" << n_distinct_refs << " referenced nodes)..." << std::endl;

  {
    std::ifstream rc_f(refcount_path.c_str(), std::ios::binary);
    std::ofstream nc_f(coords_bin_path.c_str(), std::ios::binary);
    if (!nc_f) Rcpp::stop("Cannot create temp file: " + coords_bin_path);

    int64_t rc_id = 0; uint32_t rc_count = 0;
    auto read_rc = [&](){
      return static_cast<bool>(rc_f.read(reinterpret_cast<char*>(&rc_id), sizeof(int64_t)) &&
                                rc_f.read(reinterpret_cast<char*>(&rc_count), sizeof(uint32_t)));
    };
    bool rc_has = read_rc();

    int64_t last_node_id = -1;
    bool first_node = true;

    stream_osmpbf(path,
      [&](const OsmNode &node){
        if (!first_node && node.id <= last_node_id) {
          Rcpp::stop("cpp_extract disk mode requires nodes to appear in strictly increasing id order in the "
                     ".pbf file (standard for osmium/osmconvert/Geofabrik extracts) -- this file violates that "
                     "assumption. Omit max_memory to use the in-RAM path instead.");
        }
        last_node_id = node.id; first_node = false;

        while (rc_has && rc_id < node.id) rc_has = read_rc();

        if (rc_has && rc_id == node.id) {
          float lat = static_cast<float>(node.lat), lon = static_cast<float>(node.lon);
          nc_f.write(reinterpret_cast<const char*>(&rc_id), sizeof(int64_t));
          nc_f.write(reinterpret_cast<const char*>(&lat), sizeof(float));
          nc_f.write(reinterpret_cast<const char*>(&lon), sizeof(float));
          nc_f.write(reinterpret_cast<const char*>(&rc_count), sizeof(uint32_t));
          rc_has = read_rc();
        }
      },
      std::function<void(const OsmWay&)>());
  }
  std::remove(refcount_path.c_str());

  // ---- Pass 4/5: (node_id, way_idx, position) sorted by node_id, merge-joined against node_coords -> (way_idx,position,coord,ref_count) sorted by (way_idx,position) ----
  if (progress) Rcpp::Rcout << "pass 4/5: attaching coordinates to each way's node sequence..." << std::endl;

  ExternalSorter<IdPosRecord, CompareIdPosById> idpos_sorter(idpos_chunk, temp_dir, "cppRouting_osm_idpos");
  {
    WaysFileReader ways_reader(ways_path);
    double speed; int direction, clazz; std::vector<int64_t> refs;
    int64_t way_idx = 0;
    while (ways_reader.read_way(speed, direction, clazz, refs)){
      for (size_t p = 0; p < refs.size(); p++)
        idpos_sorter.add(IdPosRecord{refs[p], way_idx, static_cast<int32_t>(p), 0});
      way_idx++;
    }
  }
  idpos_sorter.finish();

  ExternalSorter<WayPosRecord, CompareWayPos> waypos_sorter(waypos_chunk, temp_dir, "cppRouting_osm_waypos");
  {
    std::ifstream nc_f(coords_bin_path.c_str(), std::ios::binary);
    int64_t nc_id = 0; float nc_lat = 0, nc_lon = 0; uint32_t nc_count = 0;
    auto read_nc = [&](){
      return static_cast<bool>(nc_f.read(reinterpret_cast<char*>(&nc_id), sizeof(int64_t)) &&
                                nc_f.read(reinterpret_cast<char*>(&nc_lat), sizeof(float)) &&
                                nc_f.read(reinterpret_cast<char*>(&nc_lon), sizeof(float)) &&
                                nc_f.read(reinterpret_cast<char*>(&nc_count), sizeof(uint32_t)));
    };
    bool nc_has = read_nc();

    for (auto it = idpos_sorter.iterator(); it.has_next(); ){
      IdPosRecord rec = it.next();
      while (nc_has && nc_id < rec.node_id) nc_has = read_nc();
      if (nc_has && nc_id == rec.node_id)
        waypos_sorter.add(WayPosRecord{rec.way_idx, rec.position, nc_lat, nc_lon, nc_count});
      // else: node has no coordinate (extract boundary) -- omitted; the
      // final grouping pass detects the resulting short group and drops
      // the whole way, same policy as the in-RAM pipeline.
    }
  }
  waypos_sorter.finish();
  std::remove(coords_bin_path.c_str());

  // ---- Pass 5/5: group by way_idx (already contiguous, sorted by (way_idx,position)), split at junctions, write edges.csv; collect graph-node coords for dedup ----
  if (progress) Rcpp::Rcout << "pass 5/5: splitting ways at junctions and writing edges..." << std::endl;

  std::ofstream edges_f(edges_path.c_str());
  if (!edges_f) Rcpp::stop("Cannot open for writing: " + edges_path);
  edges_f << "from,to,dist\n" << std::setprecision(10);

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

  std::unique_ptr<ExternalSorter<CoordRecord, CompareCoordById> > coord_sorter;
  if (want_coords) coord_sorter.reset(new ExternalSorter<CoordRecord, CompareCoordById>(coord_chunk, temp_dir, "cppRouting_osm_coordout"));

  long long n_edges = 0;
  int n_boundary_ways_skipped = 0;
  {
    WaysFileReader ways_reader(ways_path);
    auto wp_it = waypos_sorter.iterator();
    bool wp_pending = wp_it.has_next();
    WayPosRecord wp_next = wp_pending ? wp_it.next() : WayPosRecord{};

    std::vector<WayPosRecord> group;
    double speed; int direction, clazz; std::vector<int64_t> refs;
    int64_t way_idx = 0;

    while (ways_reader.read_way(speed, direction, clazz, refs)){
      group.clear();
      while (wp_pending && wp_next.way_idx == way_idx){
        group.push_back(wp_next);
        wp_pending = wp_it.has_next();
        if (wp_pending) wp_next = wp_it.next();
      }

      if (group.size() == refs.size()){
        size_t n = group.size();
        size_t seg_start = 0;
        for (size_t i = 1; i < n; i++){
          bool is_last = (i == n - 1);
          bool is_junction = group[i].ref_count >= 2;
          if (!is_last && !is_junction) continue;

          int64_t from_osm = refs[seg_start];
          int64_t to_osm = refs[i];
          double length_m = 0;
          for (size_t j = seg_start; j < i; j++)
            length_m += osm_tags::haversine_m(group[j].lat, group[j].lon, group[j + 1].lat, group[j + 1].lon);

          if (from_osm != to_osm){
            double length_km = length_m / 1000.0;
            double cost = 60.0 * length_km / speed;
            if (want_coords) {
              coord_sorter->add(CoordRecord{from_osm, group[seg_start].lat, group[seg_start].lon});
              coord_sorter->add(CoordRecord{to_osm, group[i].lat, group[i].lon});
            }
            auto write_aux_row = [&](){
              if (!want_aux) return;
              if (keep_clazz && keep_length) aux_f << clazz << "," << length_km << "\n";
              else if (keep_clazz) aux_f << clazz << "\n";
              else aux_f << length_km << "\n";
            };
            // ".0" forces from/to to read back as plain double, not
            // bit64::integer64 (see write_edges_from_ways_csv() in
            // io_osm_graph.cpp for why: real OSM ids exceed 2^31, and
            // integer64 without the bit64 package holds outright garbage,
            // not just a display quirk).
            if (direction != 2){ edges_f << from_osm << ".0," << to_osm << ".0," << cost << "\n"; write_aux_row(); n_edges++; }
            if (direction != 1){ edges_f << to_osm << ".0," << from_osm << ".0," << cost << "\n"; write_aux_row(); n_edges++; }
          }

          seg_start = i;
        }
      } else {
        n_boundary_ways_skipped++;
      }

      way_idx++;
    }
  }
  std::remove(ways_path.c_str());

  long long n_nodes = 0;
  if (want_coords) {
    coord_sorter->finish();

    // Dedup the graph-node coordinates (each node was added once per way
    // that uses it as a split-point endpoint, so duplicates across ways are
    // expected) and write the final coords.csv.
    std::ofstream coords_f(coords_path.c_str());
    if (!coords_f) Rcpp::stop("Cannot open for writing: " + coords_path);
    coords_f << "node_id,X,Y\n" << std::setprecision(9);

    auto cit = coord_sorter->iterator();
    if (cit.has_next()){
      CoordRecord cur = cit.next();
      while (true){
        bool has_next = cit.has_next();
        CoordRecord nxt = has_next ? cit.next() : CoordRecord{};
        if (!has_next || nxt.node_id != cur.node_id){
          coords_f << cur.node_id << ".0," << cur.lon << "," << cur.lat << "\n"; // see the ".0" comment in the edges_f writer above
          n_nodes++;
          if (!has_next) break;
        }
        cur = nxt;
      }
    }
  }

  if (progress) {
    if (want_coords) {
      Rcpp::Rcout << "done: " << n_nodes << " nodes, " << n_edges << " edges ("
                  << n_boundary_ways_skipped << " way(s) skipped, crossing the extract's boundary) from "
                  << n_kept_ways << " kept ways" << std::endl;
    } else {
      Rcpp::Rcout << "done: " << n_edges << " edges (coordinates not extracted) ("
                  << n_boundary_ways_skipped << " way(s) skipped, crossing the extract's boundary) from "
                  << n_kept_ways << " kept ways" << std::endl;
    }
  }

  return Rcpp::List::create(Rcpp::Named("nbnode") = want_coords ? static_cast<double>(n_nodes) : NA_REAL,
                             Rcpp::Named("nbedge") = static_cast<double>(n_edges));
}

// [[Rcpp::export]]
Rcpp::List cpp_build_osm_graph_disk_csv(std::string path,
                                          std::vector<std::string> profile_highway,
                                          std::vector<int> profile_clazz,
                                          std::vector<double> profile_speed,
                                          bool use_maxspeed,
                                          bool assume_oneway_motorway,
                                          std::string mode,
                                          double max_memory_mb,
                                          std::string temp_dir,
                                          std::string edges_path,
                                          std::string coords_path,
                                          std::string aux_path,
                                          bool want_coords,
                                          bool keep_clazz,
                                          bool keep_length,
                                          bool progress){
  return build_osm_graph_disk(path, profile_highway, profile_clazz, profile_speed,
                               use_maxspeed, assume_oneway_motorway, mode, max_memory_mb, temp_dir,
                               edges_path, coords_path, aux_path, want_coords, keep_clazz, keep_length, progress);
}
