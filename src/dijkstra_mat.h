#ifndef DIJKSTRA_MAT_H_INCLUDED
#define DIJKSTRA_MAT_H_INCLUDED

#include <Rcpp.h>

Rcpp::NumericMatrix InternalDijkstra_mat(std::vector<int> gfrom,std::vector<int> gto,std::vector<double> gw,int NbNodes,std::vector<int> dep, std::vector<int> arr);

#endif
