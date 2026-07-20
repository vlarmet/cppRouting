#include "external_sort.h"
#include <Rcpp.h>
#include <algorithm>

// Test-only export: sorts a vector of doubles (as int64_t) through the
// external sorter with an artificially small chunk size (to force multiple
// runs + a real k-way merge even on tiny inputs) and returns the result, so
// the R test suite can compare it against sort(). Not part of the public
// package surface.
// [[Rcpp::export]]
Rcpp::NumericVector cpp_test_external_sort(std::vector<double> values, int chunk_capacity, std::string temp_dir){
  ExternalSorter<int64_t> sorter(static_cast<size_t>(chunk_capacity), temp_dir, "test_extsort");
  for (double v : values) sorter.add(static_cast<int64_t>(v));
  sorter.finish();

  std::vector<double> out;
  out.reserve(values.size());
  for (auto it = sorter.iterator(); it.has_next(); ) out.push_back(static_cast<double>(it.next()));

  return Rcpp::wrap(out);
}
