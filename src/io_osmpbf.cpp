#include "io_osmpbf.h"
#include <Rcpp.h>
#include <zlib.h>
#include <fstream>
#include <cstring>
#include <stdexcept>

// Minimal hand-rolled protobuf reader, scoped to exactly what the OSM PBF
// schema (fileformat.proto / osmformat.proto) needs: varints, zigzag
// sint32/64, length-delimited bytes/embedded messages, and packed repeated
// scalar fields. No general .proto support (no maps/oneofs/groups needed).
namespace {

class PBReader {
public:
  PBReader(const unsigned char *data, size_t len) : d(data), n(len), p(0) {}

  bool eof() const { return p >= n; }

  uint64_t varint(){
    uint64_t result = 0;
    int shift = 0;
    while (true){
      if (p >= n) throw std::runtime_error("truncated varint in PBF stream");
      unsigned char b = d[p++];
      result |= static_cast<uint64_t>(b & 0x7F) << shift;
      if (!(b & 0x80)) break;
      shift += 7;
    }
    return result;
  }

  void tag(int &field_no, int &wire_type){
    uint64_t t = varint();
    field_no = static_cast<int>(t >> 3);
    wire_type = static_cast<int>(t & 0x7);
  }

  std::pair<const unsigned char*, size_t> bytes(){
    uint64_t len = varint();
    if (p + len > n) throw std::runtime_error("truncated length-delimited PBF field");
    const unsigned char *start = d + p;
    p += static_cast<size_t>(len);
    return std::make_pair(start, static_cast<size_t>(len));
  }

  void skip(int wire_type){
    switch (wire_type){
      case 0: varint(); break;
      case 1: p += 8; break;
      case 2: { uint64_t len = varint(); p += static_cast<size_t>(len); break; }
      case 5: p += 4; break;
      default: throw std::runtime_error("unsupported PBF wire type");
    }
  }

private:
  const unsigned char *d;
  size_t n;
  size_t p;
};

inline int64_t zigzag(uint64_t v){
  return static_cast<int64_t>(v >> 1) ^ -static_cast<int64_t>(v & 1);
}

std::vector<int64_t> read_packed_svarint(const unsigned char *d, size_t n){
  PBReader r(d, n);
  std::vector<int64_t> out;
  while (!r.eof()) out.push_back(zigzag(r.varint()));
  return out;
}

std::vector<uint64_t> read_packed_varint(const unsigned char *d, size_t n){
  PBReader r(d, n);
  std::vector<uint64_t> out;
  while (!r.eof()) out.push_back(r.varint());
  return out;
}

std::vector<std::string> parse_stringtable(const unsigned char *d, size_t n){
  PBReader r(d, n);
  std::vector<std::string> out;
  while (!r.eof()){
    int fn, wt; r.tag(fn, wt);
    if (fn == 1 && wt == 2){
      auto b = r.bytes();
      out.emplace_back(reinterpret_cast<const char*>(b.first), b.second);
    } else r.skip(wt);
  }
  return out;
}

void parse_dense_nodes(const unsigned char *d, size_t n, int32_t granularity,
                        int64_t lat_offset, int64_t lon_offset,
                        const std::function<void(const OsmNode&)> &on_node){
  if (!on_node) return;

  PBReader r(d, n);
  std::vector<int64_t> ids, lats, lons;
  while (!r.eof()){
    int fn, wt; r.tag(fn, wt);
    if (fn == 1 && wt == 2){ auto b = r.bytes(); ids = read_packed_svarint(b.first, b.second); }
    else if (fn == 8 && wt == 2){ auto b = r.bytes(); lats = read_packed_svarint(b.first, b.second); }
    else if (fn == 9 && wt == 2){ auto b = r.bytes(); lons = read_packed_svarint(b.first, b.second); }
    else r.skip(wt); // denseinfo(5), keys_vals(10) -- node tags not needed for the road graph
  }

  int64_t id = 0, lat = 0, lon = 0;
  size_t count = ids.size();
  for (size_t i = 0; i < count; i++){
    id += ids[i];
    lat += (i < lats.size() ? lats[i] : 0);
    lon += (i < lons.size() ? lons[i] : 0);
    OsmNode node;
    node.id = id;
    node.lat = 1e-9 * (static_cast<double>(lat_offset) + static_cast<double>(granularity) * lat);
    node.lon = 1e-9 * (static_cast<double>(lon_offset) + static_cast<double>(granularity) * lon);
    on_node(node);
  }
}

void parse_plain_node(const unsigned char *d, size_t n, int32_t granularity,
                       int64_t lat_offset, int64_t lon_offset,
                       const std::function<void(const OsmNode&)> &on_node){
  if (!on_node) return;

  PBReader r(d, n);
  OsmNode node; node.id = 0; node.lat = 0; node.lon = 0;
  int64_t raw_lat = 0, raw_lon = 0;
  while (!r.eof()){
    int fn, wt; r.tag(fn, wt);
    if (fn == 1 && wt == 0) node.id = zigzag(r.varint());
    else if (fn == 8 && wt == 0) raw_lat = zigzag(r.varint());
    else if (fn == 9 && wt == 0) raw_lon = zigzag(r.varint());
    else r.skip(wt);
  }
  node.lat = 1e-9 * (static_cast<double>(lat_offset) + static_cast<double>(granularity) * raw_lat);
  node.lon = 1e-9 * (static_cast<double>(lon_offset) + static_cast<double>(granularity) * raw_lon);
  on_node(node);
}

void parse_way(const unsigned char *d, size_t n, const std::vector<std::string> &stringtable,
               const std::function<void(const OsmWay&)> &on_way){
  if (!on_way) return;

  PBReader r(d, n);
  OsmWay way; way.id = 0;
  std::vector<uint64_t> keys, vals;
  while (!r.eof()){
    int fn, wt; r.tag(fn, wt);
    if (fn == 1 && wt == 0) way.id = static_cast<int64_t>(r.varint());
    else if (fn == 2 && wt == 2){ auto b = r.bytes(); keys = read_packed_varint(b.first, b.second); }
    else if (fn == 3 && wt == 2){ auto b = r.bytes(); vals = read_packed_varint(b.first, b.second); }
    else if (fn == 8 && wt == 2){
      auto b = r.bytes();
      auto deltas = read_packed_svarint(b.first, b.second);
      int64_t ref = 0;
      way.refs.reserve(deltas.size());
      for (int64_t dref : deltas){ ref += dref; way.refs.push_back(ref); }
    }
    else r.skip(wt); // info(4) -- metadata not needed
  }

  size_t n_tags = std::min(keys.size(), vals.size());
  way.tags.reserve(n_tags);
  for (size_t i = 0; i < n_tags; i++){
    if (keys[i] < stringtable.size() && vals[i] < stringtable.size())
      way.tags.push_back(OsmTag{stringtable[keys[i]], stringtable[vals[i]]});
  }
  on_way(way);
}

void parse_primitive_group(const unsigned char *d, size_t n, const std::vector<std::string> &stringtable,
                            int32_t granularity, int64_t lat_offset, int64_t lon_offset,
                            const std::function<void(const OsmNode&)> &on_node,
                            const std::function<void(const OsmWay&)> &on_way){
  PBReader r(d, n);
  while (!r.eof()){
    int fn, wt; r.tag(fn, wt);
    if (fn == 1 && wt == 2){ auto b = r.bytes(); parse_plain_node(b.first, b.second, granularity, lat_offset, lon_offset, on_node); }
    else if (fn == 2 && wt == 2){ auto b = r.bytes(); parse_dense_nodes(b.first, b.second, granularity, lat_offset, lon_offset, on_node); }
    else if (fn == 3 && wt == 2){ auto b = r.bytes(); parse_way(b.first, b.second, stringtable, on_way); }
    else r.skip(wt); // relations(4), changesets(5) -- turn restrictions etc., out of scope for now
  }
}

void parse_primitive_block(const unsigned char *d, size_t n,
                            const std::function<void(const OsmNode&)> &on_node,
                            const std::function<void(const OsmWay&)> &on_way){
  PBReader r(d, n);
  std::vector<std::string> stringtable;
  int32_t granularity = 100;
  int64_t lat_offset = 0, lon_offset = 0;
  std::vector<std::pair<const unsigned char*, size_t> > groups;

  while (!r.eof()){
    int fn, wt; r.tag(fn, wt);
    if (fn == 1 && wt == 2){ auto b = r.bytes(); stringtable = parse_stringtable(b.first, b.second); }
    else if (fn == 2 && wt == 2){ groups.push_back(r.bytes()); }
    else if (fn == 17 && wt == 0){ granularity = static_cast<int32_t>(r.varint()); }
    else if (fn == 19 && wt == 0){ lat_offset = static_cast<int64_t>(r.varint()); }
    else if (fn == 20 && wt == 0){ lon_offset = static_cast<int64_t>(r.varint()); }
    else r.skip(wt);
  }

  for (auto &g : groups)
    parse_primitive_group(g.first, g.second, stringtable, granularity, lat_offset, lon_offset, on_node, on_way);
}

struct BlobHeaderInfo { std::string type; int32_t datasize; };

BlobHeaderInfo parse_blob_header(const unsigned char *d, size_t n){
  PBReader r(d, n);
  BlobHeaderInfo info; info.datasize = 0;
  while (!r.eof()){
    int fn, wt; r.tag(fn, wt);
    if (fn == 1 && wt == 2){ auto b = r.bytes(); info.type.assign(reinterpret_cast<const char*>(b.first), b.second); }
    else if (fn == 3 && wt == 0){ info.datasize = static_cast<int32_t>(r.varint()); }
    else r.skip(wt);
  }
  return info;
}

std::vector<unsigned char> parse_blob(const unsigned char *d, size_t n){
  PBReader r(d, n);
  const unsigned char *raw_ptr = nullptr; size_t raw_len = 0;
  const unsigned char *zlib_ptr = nullptr; size_t zlib_len = 0;
  int32_t raw_size = 0;

  while (!r.eof()){
    int fn, wt; r.tag(fn, wt);
    if (fn == 1 && wt == 2){ auto b = r.bytes(); raw_ptr = b.first; raw_len = b.second; }
    else if (fn == 2 && wt == 0){ raw_size = static_cast<int32_t>(r.varint()); }
    else if (fn == 3 && wt == 2){ auto b = r.bytes(); zlib_ptr = b.first; zlib_len = b.second; }
    else r.skip(wt); // lzma(4)/obsolete bzip2(5)/lz4(6)/zstd(7) -- not produced by standard extracts
  }

  if (zlib_ptr){
    std::vector<unsigned char> out(raw_size > 0 ? static_cast<size_t>(raw_size) : zlib_len * 4);
    uLongf destlen = static_cast<uLongf>(out.size());
    int ret = uncompress(out.data(), &destlen, zlib_ptr, static_cast<uLong>(zlib_len));
    if (ret != Z_OK) Rcpp::stop("zlib decompression failed for a PBF blob (corrupt file?)");
    out.resize(destlen);
    return out;
  }
  if (raw_ptr) return std::vector<unsigned char>(raw_ptr, raw_ptr + raw_len);
  Rcpp::stop("PBF blob has no raw or zlib payload (lzma/bzip2/lz4/zstd-compressed PBF files are not supported)");
}

}

void stream_osmpbf(const std::string &path,
                    const std::function<void(const OsmNode&)> &on_node,
                    const std::function<void(const OsmWay&)> &on_way){
  std::ifstream f(path.c_str(), std::ios::binary);
  if (!f) Rcpp::stop("Cannot open file: " + path);

  while (true){
    unsigned char len_buf[4];
    f.read(reinterpret_cast<char*>(len_buf), 4);
    std::streamsize got = f.gcount();
    if (got == 0) break;
    if (got != 4) Rcpp::stop("Truncated PBF blob header length in: " + path);

    uint32_t header_len = (static_cast<uint32_t>(len_buf[0]) << 24) | (static_cast<uint32_t>(len_buf[1]) << 16) |
                           (static_cast<uint32_t>(len_buf[2]) << 8)  |  static_cast<uint32_t>(len_buf[3]);

    std::vector<unsigned char> header_buf(header_len);
    f.read(reinterpret_cast<char*>(header_buf.data()), header_len);
    if (static_cast<uint32_t>(f.gcount()) != header_len) Rcpp::stop("Truncated PBF blob header in: " + path);

    BlobHeaderInfo bh = parse_blob_header(header_buf.data(), header_len);

    std::vector<unsigned char> blob_buf(static_cast<size_t>(bh.datasize));
    f.read(reinterpret_cast<char*>(blob_buf.data()), bh.datasize);
    if (f.gcount() != bh.datasize) Rcpp::stop("Truncated PBF blob body in: " + path);

    if (bh.type == "OSMData"){
      std::vector<unsigned char> block = parse_blob(blob_buf.data(), blob_buf.size());
      parse_primitive_block(block.data(), block.size(), on_node, on_way);
    }
    // "OSMHeader" blob (always first): bbox/writer metadata, not needed here.
  }
}

// [[Rcpp::export]]
Rcpp::List cpp_osmpbf_counts(std::string path){
  double n_nodes = 0, n_ways = 0, n_tagged_ways = 0;
  double lat_sum = 0, lon_sum = 0;

  stream_osmpbf(path,
    [&](const OsmNode &node){ n_nodes++; lat_sum += node.lat; lon_sum += node.lon; },
    [&](const OsmWay &way){ n_ways++; if (!way.tags.empty()) n_tagged_ways++; });

  return Rcpp::List::create(
    Rcpp::Named("n_nodes") = n_nodes,
    Rcpp::Named("n_ways") = n_ways,
    Rcpp::Named("n_tagged_ways") = n_tagged_ways,
    Rcpp::Named("mean_lat") = lat_sum / n_nodes,
    Rcpp::Named("mean_lon") = lon_sum / n_nodes
  );
}
