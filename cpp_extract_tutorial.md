# Building routable graphs with `cpp_extract()`

`cpp_extract()` is a single entry point that turns a geographic linestring source into a
[`makegraph()`](R/graphs.R)-ready graph, without needing GDAL.

Two families of input are supported:
- **OpenStreetMap**: `.osm.pbf` files, an `osmdata` object, an `osmextract`-cached extract
- **generic `LINESTRING`/`MULTILINESTRING` dataset**: an `sf` object, or a shapefile.
  **Assumed to already be topologically clean.**

It scales from a handful of roads up to a full country or continent without external
dependencies.

Five kinds of input are recognized automatically:

| Input | How you pass it |
|---|---|
| Shapefile of linestrings (e.g. IGN ROUTE500, BD TOPO) | a path to the `.shp` file |
| In-memory `sf` object (`LINESTRING`/`MULTILINESTRING`) | the `sf` object itself |
| OpenStreetMap extract | a path to a `.osm.pbf` file |
| `osmdata` query result | the full result of `osmdata::osmdata_sf()` |
| `osmextract`-cached extract | the raw cached `.pbf` path (see below) |

GeoPackage is **not** supported directly -- convert it to a shapefile first (`sf::st_write(sf::st_read("file.gpkg"), "file.shp")`) and pass that instead.

This tutorial builds a tiny synthetic road network so every example below is copy-pasteable
and runs in a few milliseconds, then shows the real-world patterns (shapefile on disk, OSM PBF,
`osmdata`, `osmextract`) and what changes at country/continent scale.

```r
library(cppRouting)
library(sf)
```

## How the arguments are organized

`cpp_extract()` has a lot of arguments, but most of them only apply to *one* of the two input
families -- shapefile/`sf` on one side, OSM (`.pbf`/`osmdata`/`osmextract`) on the other. The
function signature groups them the same way: family-specific arguments are elements of a
`geom_params = list(...)` or `osm_params = list(...)` argument, rather than individual top-level arguments.
This tutorial follows the same grouping, so you only need to read the part that matches your
source:

| | `geom_params = list(...)` (shapefile / `sf`) | `osm_params = list(...)` (`.pbf`/`osmdata`/`osmextract`) |
|---|---|---|
| **Cost** | `cost_field` OR `speed_field` (+ `length_field`) | derived from `profile` + `use_maxspeed` |
| **Direction** | `direction_field` + `direction_values`, or `directed` | derived from OSM's `oneway` tag + `assume_oneway_motorway` |
| **Extra attributes** | `keep_fields` | `keep_clazz` |
| **Scale** | -- | `max_mem`, `tmp_path` (`.pbf` only) |

**Shared by every source, regardless of family** -- these stay as individual top-level
arguments, not nested in `osm_params`/`geom_params`: `output`/`out_prefix`, `coords`/`epsg`, `keep_length`,
and `snap_digits` (endpoint dedup for shapefile/`sf`; matches vertices against `$osm_points` for
`osmdata`; unused for `.pbf`, which has stable node ids -- the one argument that doesn't fit
cleanly in either bucket, so it's left flat). These are covered first, once, below -- rather
than repeated in both the shapefile/`sf` and OSM sections.

## Arguments shared by every source

### Output modes: `output = "graph" | "data.frame" | "csv"`

- **`"graph"`** (default): returns a ready-to-use `cppRouting_graph` (internally calls
  `makegraph()` for you).
- **`"data.frame"`**: returns `list(data, coords, aux)` -- the edges/coords/attribute tables,
  without building the graph. Useful to inspect or post-process before calling `makegraph()`
  yourself.
- **`"csv"`** (needs `out_prefix`): writes `<out_prefix>_edges.csv` and, if `coords = TRUE`,
  `<out_prefix>_coords.csv`, and returns their paths invisibly. For OSM sources in particular,
  this streams edges/coordinates straight to disk as they're produced instead of building the
  full result in R memory first -- see [Scaling up](#scaling-up-country-and-continent-size-extracts).
  `csv` always carries the same information as `data.frame`: if `keep_fields`/`keep_clazz`/
  `keep_length` produces an `aux` table, it's written to `<out_prefix>_aux.csv` too.

```r
res <- cpp_extract(roads, geom_params = list(speed_field = "vitesse"), output = "data.frame")
res$data    #  from, to, dist
res$coords  #  node_id, X, Y

paths <- cpp_extract(roads, geom_params = list(speed_field = "vitesse"), output = "csv", out_prefix = "network")
#> paths[1] == "network_edges.csv", paths[2] == "network_coords.csv"
```

### Node coordinates: `coords`, `epsg`

By default `cpp_extract()` also extracts each node's coordinates (`coords = TRUE`), attached
to the graph as `g$coords` (or `res$coords` / `<out_prefix>_coords.csv`, depending on `output`)
and usable by `makegraph()`'s A\*/NBA algorithms.

**Geographic (lat/lon) coordinates need `epsg`.** `get_distance_pair()`/`get_path_pair()`'s
A\*/NBA algorithms use euclidean distance between node coordinates as a heuristic, which only
approximates the real edge cost if those coordinates are in a *projected* (meter-scale) CRS.
OSM sources (`.pbf`/`osmdata`) are always geographic; a shapefile/`sf` source is geographic if
its `.prj`/CRS says so. If `coords = TRUE` (the default) and the source is geographic,
`cpp_extract()` refuses outright rather than silently handing back coordinates that would
break those algorithms:

```r
cpp_extract("germany-latest.osm.pbf")
#> Error: Coordinates are geographic (lat/lon degrees), not projected. [...]
#> Either set coords = FALSE if you don't need node coordinates, or supply
#> epsg = <a projected EPSG code> to reproject them (requires the 'sf' package).
```

Two ways out:

```r
# Option 1: you don't need coordinates at all (e.g. Dijkstra/bidirectional only,
# or you already have coordinates from elsewhere)
g <- cpp_extract("germany-latest.osm.pbf", coords = FALSE)

# Option 2: reproject to a projected CRS (requires the 'sf' package)
g <- cpp_extract("germany-latest.osm.pbf", epsg = 25832)  # UTM zone 32N, e.g.
```

`epsg` accepts any EPSG code (or anything `sf::st_crs()` understands) and is silently ignored
if the source is already projected. It is *not* supported together with `max_mem` (OSM-only,
see [Scaling up](#scaling-up-country-and-continent-size-extracts)) -- reprojecting there would
mean reading the whole coords table back into R, defeating the point of that pipeline's bounded
RAM.

### Time *and* distance: `keep_length`

A common need: route by time, but still report the physical distance of the chosen path.
`keep_length = TRUE` keeps each edge's length (in km) as `aux`, computed independently of how
the main `dist` (cost) was resolved -- for shapefile/`sf` sources this is the geometric length,
for OSM sources the haversine length:

```r
res <- cpp_extract(roads,
                    geom_params = list(speed_field      = "vitesse",
                                   direction_field  = "sens",
                                   direction_values = list(twoway = "Double sens", forward = "Sens unique")),
                    keep_length = TRUE, output = "data.frame")
cbind(res$data, res$aux)
#>   from to     dist   length
#> 1    0  1 1.200000 1.000000
#> 2    1  0 1.200000 1.000000
#> 3    1  2 1.200000 1.000000
#> 4    2  0 2.828427 1.414214
#> 5    0  2 2.828427 1.414214
```

For shapefile/`sf` sources with a geographic CRS, this is the raw degree-based length
(unreliable as a physical distance either way -- reproject with `epsg`, or use `keep_fields`
with your own length column instead).

For OSM sources, `keep_fields` doesn't apply (cost/direction come from tags, not attribute
columns you name -- see the table above), so `keep_length` is the only way to get physical
distance alongside a time-based OSM graph:

```r
g <- cpp_extract("bourgogne.osm.pbf", epsg = 2154, keep_length = TRUE, output = "data.frame")
g$aux$length  # km, one row per edge in g$data
```

Like `keep_fields`/`keep_clazz`, `keep_length` is written to `<out_prefix>_aux.csv` when
`output = "csv"` -- including through the `max_mem` disk-backed pipeline, where `clazz`/
`length` are already computed at exactly the point each edge is written, so this doesn't cost
an extra pass (unlike `coords`, which does need one).

## Non-OSM sources: shapefile and `sf`

For a shapefile or `sf` source, `cpp_extract()` assumes the network is already **topologically
clean**: lines are expected to only meet at their endpoints, and topology is reconstructed purely
from coordinates -- there is no stable node id to rely on (unlike OSM). Concretely, that means
three edges like this:

```
 (0,0)                      (200,0)
   •───────────────────────────•
             Edge 1            │
                               │
                               │ Edge 2
                               │
                               •───────────────────────────•
                            (200,150)          Edge 3      (400,150)
```
*Figure 1 -- three edges forming a chain. At a glance, Edge 2 and Edge 3 look like they share the
same point at `(200,150)`.*

In real digitized data, "the same point" is rarely bit-for-bit identical -- two features drawn
(or exported) independently can end up a fraction of a millimeter apart even though they
represent the same real-world junction. Zooming in on that Edge 2 / Edge 3 junction:

```
        Edge 2 ends here
              │
              ▼
              •   (200.00000, 150.00000)
              ┆
              ┆   0.00002 units apart (20 micrometers, invisible at Figure 1's scale)
              ┆
              •   (200.00002, 150.00000)
              ▲
              │
        Edge 3 starts here
```
*Figure 2 -- zoomed in on the junction: the two endpoints are not perfectly coincident.*

This is exactly what `snap_digits` is for: coordinates are rounded to `snap_digits` decimal
places before deduplication, so two endpoints that round to the same value become **one** node.
With the default `snap_digits = 4` for projected (meter-scale) coordinates, both endpoints above
round to `(200.0000, 150.0000)` and correctly collapse into a single node -- too fine (too many
decimal places) a `snap_digits` (e.g. `7`) would keep them distinct, silently splitting one real
junction into two disconnected nodes:

```
 node 0                node 1                 node 2                node 3
 (0,0)                 (200,0)                (200,150)              (400,150)
   •──────────────────────•──────────────────────•──────────────────────•
        Edge 1: (0,1)            Edge 2: (1,2)          Edge 3: (2,3)
```
*Figure 3 -- the extracted nodes: 3 edges, 4 nodes. Node 2 is exactly the Edge 2/Edge 3 junction
from Figure 2 -- thanks to snapping, it's one node (not two), so Edge 2 and Edge 3 connect
properly. Node ids are assigned in encounter order (`0, 1, 2, ...`), ready for `from`/`to` in the
edge table.*

Left at its default (`NULL`), `snap_digits` resolves to `7` for geographic coordinates or `4`
for projected ones -- geographic coordinates need many decimal places to distinguish nearby
points, but that same precision is far too tight for projected (meter-scale) coordinates, where
it starts picking up floating-point noise between independently-digitized but
logically-identical points rather than real geometric differences (validated against IGN's own
ROUTE500 node layer: exact-match rates collapsed from 99.9% at `snap_digits = 4` to 0% at
`snap_digits = 6`, purely from sub-micrometer noise). Lower it further if nearby-but-distinct
endpoints from imprecise digitizing are wrongly being merged into one node.

### A minimal example

Three streets forming a triangle, with a direction column (`sens`, French for "direction") and
a speed column (`vitesse`, km/h), in a projected CRS (Lambert-93, the standard French projection --
any projected, meter-based CRS works the same way):

```r
roads <- st_sf(
  name    = c("Rue A", "Rue B", "Rue C"),
  sens    = c("Double sens", "Sens unique", "Double sens"),
  vitesse = c(50, 50, 30),
  geometry = st_sfc(
    st_linestring(matrix(c(0,0, 1000,0),       ncol = 2, byrow = TRUE)),  # Rue A: 1000m
    st_linestring(matrix(c(1000,0, 1000,1000), ncol = 2, byrow = TRUE)),  # Rue B: 1000m
    st_linestring(matrix(c(1000,1000, 0,0),    ncol = 2, byrow = TRUE)),  # Rue C: ~1414m
    crs = 2154
  )
)

g <- cpp_extract(roads,
                  geom_params = list(speed_field      = "vitesse",
                                 direction_field  = "sens",
                                 direction_values = list(twoway = "Double sens", forward = "Sens unique")))

g$dict
#>   ref id
#> 1   0  0
#> 2   1  1
#> 3   2  2

get_distance_pair(g, from = "0", to = "1")
#> [1] 1.2   # minutes: 1000m at 50km/h
```

`g` is a normal `cppRouting_graph` -- everything else in the package (`get_distance_pair()`,
`get_path_pair()`, `cpp_contract()`, isochrones, traffic assignment...) works on it exactly as
if you had built it by hand with `makegraph()`.

Passing a path to a `.shp` file works identically -- `cpp_extract()` dispatches on the type/
extension of `x` automatically:

```r
st_write(roads, "roads.shp")
g <- cpp_extract("roads.shp",
                  geom_params = list(speed_field      = "vitesse",
                                 direction_field  = "sens",
                                 direction_values = list(twoway = "Double sens", forward = "Sens unique")))
```

### Cost: `cost_field`, `speed_field`, `length_field`

Exactly one of these controls edge cost:

- **`cost_field`**: a column already holding a ready-to-use cost (e.g. a precomputed travel
  time). Used as-is, nothing else needed.
- **`speed_field`** (km/h) **+ `length_field`** (km): cost is computed as `60 * length / speed`
  (minutes). This mirrors IGN ROUTE500's own convention (`LONGUEUR` in km, paired with a
  vitesse column).
- **`speed_field` alone** (no `length_field`): length is instead computed directly from the
  line geometry. cppRouting converts it from the CRS's native unit (assumed meters, true for
  any typical projected CRS) to km before applying the same formula -- so `vitesse = 50` on a
  1000m line still gives `1.2` minutes, not `1200`.
- **Neither**: cost falls back to the raw geometric length (in the CRS's native unit -- meters
  for a projected CRS), i.e. a physical-distance graph rather than a time-based one.

```r
# cost_field: use a precomputed column directly
cpp_extract(roads, geom_params = list(cost_field = "vitesse"))  # (contrived: normally a time/cost column)

# speed_field + explicit length_field (km)
roads$longueur_km <- c(1, 1, 1.4142)
cpp_extract(roads, geom_params = list(speed_field = "vitesse", length_field = "longueur_km"))

# speed_field alone: length computed from geometry, unit-converted automatically
cpp_extract(roads, geom_params = list(speed_field = "vitesse"))
```

### Direction: `direction_field` / `direction_values` / `directed`

If your data has a direction column (like ROUTE500's `SENS`), map its values with
`direction_values`:

```r
direction_values = list(
  forward  = "Sens unique",    # only first point -> last point
  backward = "Sens inverse",   # only last point -> first point
  twoway   = "Double sens"     # both directions
)
```

Any value of `direction_field` not listed anywhere is treated as `twoway`, with a warning --
useful for catching a typo or an unexpected category in real data.

Without a `direction_field`, use the simpler `directed` flag: `TRUE` keeps every line as
digitized (first point to last point only), `FALSE` (the default) makes every line
bidirectional -- the safer default when you have no direction information at all.

### Extra attributes: `keep_fields`

Keep additional per-edge columns as `aux` (see `makegraph()`), e.g. a physical distance
alongside a time-based cost:

```r
roads$distance_m <- c(1000, 1000, 1414.2)
res <- cpp_extract(roads,
                    geom_params = list(speed_field      = "vitesse",
                                   direction_field  = "sens",
                                   direction_values = list(twoway = "Double sens", forward = "Sens unique"),
                                   keep_fields      = "distance_m"),
                    output = "data.frame")
res$aux
#>   distance_m
#> 1     1000.0
#> ...
```

`output = "csv"` writes the exact same `aux` to `<out_prefix>_aux.csv`, one row per edge in
`<out_prefix>_edges.csv`:

```r
paths <- cpp_extract(roads,
                      geom_params = list(speed_field      = "vitesse",
                                     direction_field  = "sens",
                                     direction_values = list(twoway = "Double sens", forward = "Sens unique"),
                                     keep_fields      = "distance_m"),
                      output = "csv", out_prefix = "network")
#> paths == c("network_edges.csv", "network_coords.csv", "network_aux.csv")
```

## OSM sources: `.pbf`, `osmdata`, `osmextract`

Unlike shapefile/`sf` sources -- which are assumed to already be topologically clean, i.e.
lines only meet at their endpoints -- OSM ways can pass through a junction in the *middle* of
a line. `cpp_extract()` detects this from OSM's own stable node ids (a node referenced by 2 or
more kept ways, or a way's own endpoint) and splits accordingly, so you don't need to
pre-process the network yourself. There is no `snap_digits`-style endpoint dedup for `.pbf`
(it doesn't need one); `osmdata` is the one exception, see below.

`from`/`to`/`node_id` for every OSM source (`.pbf`, `osmdata`, `osmextract`) are the **real OSM
node id** -- not a synthetic 0-based index -- so you can cross-reference a node directly against
`openstreetmap.org/node/<id>` or join against other OSM-derived data. For `output = "data.frame"`
they're `character`: a node id is a label, not a quantity, and this holds for **every** source
(even shapefile/`sf`, whose ids are just a synthetic index -- Module A has no native node
identity to preserve, unlike OSM) -- one consistent type for node identity, matching what
`makegraph()` itself already guarantees for `dict$ref` regardless of source. `output = "graph"`
carries the real id through exactly that `dict`: `g$dict$ref` holds it as a string, while
`g$data$from`/`g$data$to` stay `makegraph()`'s own compact integer index (its adjacency list is
always built on that, never on the original id).

`output = "csv"`: real OSM ids already exceed 2^31, so `from`/`to`/`node_id` are written with a
forced `.0` when needed (e.g. `13370063.0`) -- without it, `data.table::fread()` silently promotes
the column to `bit64::integer64` once values cross that threshold, which without the (optional)
`bit64` package installed holds garbage, not just an odd-looking print. The `.0` avoids that
ambiguity; reading it back gives a plain `double`, exact for any id up to 2^53. Note that
`fread()`'s type inference ignores quoting and leading zeros alike (a quoted `"0123"` still comes
back as the integer `123`) -- if an id's exact text ever matters, read it with
`colClasses = list(character = "node_id")`.

### Cost, direction and extra attributes: `profile`, `use_maxspeed`, `assume_oneway_motorway`, `keep_clazz`

Unlike shapefile/`sf` sources, cost and direction for OSM sources are never a column *you*
name -- they're derived automatically from OSM tags. All four below are elements of the
`osm_params = list(...)` argument, e.g. `cpp_extract(path, osm_params = list(profile = my_profile,
use_maxspeed = FALSE))`:

- **`profile`**: a `data.frame(highway, clazz, speed)` describing which `highway` tag values
  to keep and their default speed (km/h). Defaults to `osm_car_profile()`, a car profile table
  matching osm2po's own defaults. Pass your own to restrict to a subset of road types, or to
  change default speeds.
- **`use_maxspeed`** (default `TRUE`): a way's `maxspeed` tag overrides the profile's default
  speed when present and parseable (plain km/h, or `"<n> mph"`).
- **`assume_oneway_motorway`** (default `TRUE`): `motorway`/`motorway_link` ways are treated as
  one-way unless explicitly tagged `oneway=no`.
- **`keep_clazz`**: keep each edge's resolved highway class id as `aux` (the OSM equivalent of
  `keep_fields` -- see [keep_length](#time-and-distance-keep_length) for the other OSM `aux`
  option, physical distance).

A way referencing a node outside the file (a regional extract cut through a way) is dropped
entirely -- a small, expected loss limited to the extract's edge. Turn restrictions, barriers
and non-car profiles are not supported yet.

### `.osm.pbf` file

```r
g <- cpp_extract("bourgogne.osm.pbf", epsg = 2154, osm_params = list(progress = TRUE))
```

### `osmdata` result

Pass the **full** result of `osmdata::osmdata_sf()`, not just `$osm_lines`: topology
reconstruction needs `$osm_points` too (it holds the `osm_id` of every node osmdata
downloaded, including untagged vertices in the middle of a way -- exactly what lets
`cpp_extract()` recover real OSM junctions).

```r
library(osmdata)
bbox_query <- opq(bbox = "Dijon, France") |>
  add_osm_feature(key = "highway") |>
  osmdata_sf()

g <- cpp_extract(bbox_query, epsg = 2154)
```

Passing only `bbox_query$osm_lines` fails loudly (with a clear error) rather than silently
producing a graph with missing junctions.

This is also the one OSM source where `snap_digits` matters: line vertices are matched against
`$osm_points` by (rounded) coordinate, since that's the only way to recover each vertex's
`osm_id` from what `osmdata_sf()` returns. Same default (`7`, since `osmdata` results are
essentially always WGS84) and same tuning logic as the shapefile/`sf` explanation above.

### `osmextract`-cached extract

Don't let `osmextract` convert its download to `sf`/GeoPackage via GDAL first -- that
conversion only keeps *tagged* nodes by default, which loses exactly the untagged
mid-way vertices needed to reconstruct topology. Instead, get the **raw cached `.pbf`
path** and feed it to `cpp_extract()` like any other `.pbf` file:

```r
library(osmextract)
path <- oe_get("Bourgogne", download_only = TRUE, skip_vectortranslate = TRUE)
g <- cpp_extract(path, epsg = 2154)
```

**Both** `download_only = TRUE` and `skip_vectortranslate = TRUE` are required --
`download_only` alone still runs the GDAL conversion and returns a `.gpkg` path instead.

## Scaling up: country and continent-size extracts

Real numbers, measured on this machine (Ryzen-class desktop, 31GB RAM) with the in-RAM
pipeline (`output = "data.frame"`, no `max_mem`):

| Extract | Input size | Nodes | Edges | Peak RAM | Time |
|---|---|---|---|---|---|
| Bourgogne (region) | 198 MB `.osm.pbf` | 196,090 | 459,743 | ~280 MB | ~18s |
| France (country) | 5.0 GB `.osm.pbf` | 5,124,522 | 11,601,873 | ~3.5 GB | ~8 min |

`makegraph()` itself is cheaper than the extraction step in both cases (roughly half the
peak RAM) -- the extraction pipeline has to hold a much larger intermediate state (every
node *referenced* by a kept way, including non-junction pass-through nodes) before
collapsing down to the final graph's node set.

Two levers if a file doesn't comfortably fit in RAM (`.pbf` sources only):

- **`output = "csv"`** streams edges/coordinates straight to disk as they're produced,
  instead of building the full result in R memory first -- keeps peak RAM bounded by the
  first two (of three) internal passes alone. Suitable up to country scale (a few GB peak
  RAM even for all of France).
- **`max_mem`** (megabytes, requires `output = "csv"`) goes further: every intermediate
  table -- not just the final result -- is spilled to disk via external sort and
  merge-joins, keeping peak RAM close to the stated budget regardless of the input file's
  size. Meant for files that don't fit even with the `output = "csv"` optimization above
  (e.g. a full-planet extract, or a country-scale extract on a memory-constrained machine).
  Slower (more disk I/O), and requires nodes to appear in strictly increasing id order in
  the file -- standard for osmium/osmconvert/Geofabrik extracts, but violated files raise a
  clear error rather than silently producing a wrong graph.

```r
paths <- cpp_extract("europe-latest.osm.pbf",
                      osm_params = list(max_mem = 2000),  # keep peak RAM near ~2GB
                      output = "csv", out_prefix = "europe",
                      coords = FALSE)        # epsg isn't supported together with max_mem
```

If you need coordinates from a `max_mem` run, reproject `<out_prefix>_coords.csv`
yourself afterward (e.g. in chunks with `sf::st_transform()`) rather than through
`cpp_extract()` itself.

### What's multithreaded (and what isn't)

`cppRouting`'s routing algorithms (batch `get_distance_pair()`/`get_distance_matrix()` queries,
etc.) are parallelized via `RcppParallel`, controlled with `RcppParallel::setThreadOptions()`.
`cpp_extract()` itself is a different story -- checked directly against the source, only one
step is actually parallel:

| Source | Step | Multithreaded? |
|---|---|---|
| Shapefile/sf | `.shp` geometry parsing (reading/decoding each record's coordinates) | **Yes** (`RcppParallel`) |
| Shapefile | `.dbf` attribute parsing | No |
| Shapefile/sf | Coordinate dedup / edge building (after geometry is read) | No |
| `.pbf` (in-RAM or `max_mem` disk-backed) | Everything -- blob decompression, tag resolution, junction splitting, external sort/merge passes | No |
| `osmdata`/`osmextract` | Vertex matching, tag resolution, edge building | No |

In practice: more cores only speed up the `.shp` geometry-reading step of a shapefile/sf
extraction (a small fraction of total time even then, per the timings above). The entire OSM
`.pbf` pipeline -- by far the most expensive part at country/continent scale -- runs on a
single core regardless of `RcppParallel::setThreadOptions()`. This is exactly why the France
extraction above takes ~8 minutes rather than scaling down with more cores.

## Troubleshooting

- **"GeoPackage is not supported"** -- convert it first:
  `sf::st_write(sf::st_read("file.gpkg"), "file.shp")`, then pass the `.shp` path.
- **"Coordinates are geographic..."** -- see [Node coordinates](#node-coordinates-coords-epsg):
  set `coords = FALSE` or supply `epsg`.
- **"This sf object has an 'osm_id' column but isn't the full osmdata result"** -- you passed
  `bbox_query$osm_lines` (or an `osmextract`-converted `sf`/GeoPackage) instead of the full
  `osmdata` result / raw `.pbf` path. Topology reconstruction needs the point layer alongside
  the lines; treating this as a plain topologically-clean network would silently produce a
  wrong graph.
- **`download_only` alone returns a `.gpkg`, not the raw `.pbf`** -- also pass
  `skip_vectortranslate = TRUE` (see [`osmextract`-cached extract](#osmextract-cached-extract)).
- **A count of nodes/edges looks off vs. what you expected** -- for shapefile/`sf` sources,
  try adjusting `snap_digits`; for OSM sources, a way referencing a node outside the file
  (regional extract cut through a way) is dropped, which is expected at the extract's edge.

## Where to go next

Once you have a `cppRouting_graph`, everything else in the package applies as usual:

```r
get_distance_pair(g, from = "0", to = "2", algorithm = "bi")
get_path_pair(g, from = "0", to = "2")
cpp_contract(g)              # contraction hierarchies for repeated queries
get_distance_matrix(g, from = g$dict$ref, to = g$dict$ref)
```

See `?makegraph`, `?get_distance_pair` and the package's main vignette
(`vignette("cppRouting")`) for the routing side of things.
