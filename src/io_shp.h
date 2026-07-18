#ifndef IO_SHP_H
#define IO_SHP_H

#include <Rcpp.h>

// Reads the geometry part of an ESRI Shapefile (.shp), shape types PolyLine
// (3), PolyLineZ (13) and PolyLineM (23) -- Z/M values are skipped (only X/Y
// matter for graph topology). Every part of every record becomes one row
// (first point, last point, geometric length), ready to feed
// build_geom_graph() after cost/direction are resolved. part_record gives,
// for each output row, the 0-based source record index, used to replicate
// attribute columns read from the companion .dbf onto multi-part records.
Rcpp::List read_shp_geometry(std::string path);

// Reads a .dbf (dBase III, as used by shapefiles) into a data.frame. Numeric
// fields (N/F) are converted to double, everything else stays character.
Rcpp::DataFrame read_dbf(std::string path);

#endif
