#ifndef IO_OSMPBF_H
#define IO_OSMPBF_H

#include <cstdint>
#include <string>
#include <vector>
#include <functional>

struct OsmTag { std::string key; std::string value; };

struct OsmNode {
  int64_t id;
  double lat;
  double lon;
};

struct OsmWay {
  int64_t id;
  std::vector<int64_t> refs;
  std::vector<OsmTag> tags;
};

// Streams an .osm.pbf file once, invoking on_node for every node (dense or
// plain) and on_way for every way, in file order. Pass an empty
// std::function (default-constructed) for whichever callback isn't needed
// by a given pass -- the corresponding primitive group is then skipped
// without allocating/decoding its contents (only the surrounding blob is
// always decompressed, since compression applies at the blob level).
void stream_osmpbf(const std::string &path,
                    const std::function<void(const OsmNode&)> &on_node,
                    const std::function<void(const OsmWay&)> &on_way);

#endif
