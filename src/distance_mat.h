#ifndef DISTANCE_MAT_H
#define DISTANCE_MAT_H

#include <vector>
#include <RcppParallel.h>
#include <Rcpp.h>
#include "graph.h"
#include "cgraph.h"
using namespace RcppParallel;


struct distanceMat : public Worker{

  // graph pointer to avoid copying input
  const Graph *m_gr;
  IVec m_dep;
  IVec m_arr;
  bool add;

  RMatrix<double> m_result;

  // constructor
  distanceMat(Graph *gr, IVec dep, IVec arr, Rcpp::NumericMatrix result);

  // algorithms as member functions

  void dijkstra_mat(std::size_t begin, std::size_t end);


  // overload operator ()
  void operator()(std::size_t begin, std::size_t end);

};

// Same as distanceMat, generalized to n_add aux columns summed inline.
// m_result is flattened n_dep x n_arr x n_add, column-major (R array
// convention: dep varies fastest, then arr, then column) so the wrapper can
// hand it back to R as a 3-D array with dim = c(n_dep, n_arr, n_add).
struct distanceMatAdd : public Worker{

  const Graph *m_gr;
  IVec m_dep;
  IVec m_arr;
  int n_add;
  DVec m_result;

  distanceMatAdd(Graph *gr, IVec dep, IVec arr);

  void dijkstra_mat(std::size_t begin, std::size_t end);

  void operator()(std::size_t begin, std::size_t end);

};


struct distanceMatC : public Worker{

  // graph pointer to avoid copying input
  CGraph *m_gr;
  IVec m_dep;
  IVec m_arr;

  RMatrix<double> m_result;

  // constructor
  distanceMatC(CGraph *gr, IVec dep, IVec arr,
               Rcpp::NumericMatrix result);


  // overload operator ()
  void operator()(std::size_t begin, std::size_t end);

};

struct phastC : public Worker{

  // graph pointer to avoid copying input
  CGraph *m_gr;
  IVec m_dep;
  IVec m_arr;
  bool add;

  RMatrix<double> m_result;

  // constructor
  phastC(CGraph *gr, IVec dep, IVec arr,
               Rcpp::NumericMatrix result);


  // overload operator ()
  void operator()(std::size_t begin, std::size_t end);

};

// Same as phastC, generalized to n_add aux columns summed inline. m_gr->add
// / m_gr->addr are expected to already hold per-shortcut N-column sums
// (precomputed by aggCAdd), same contract as distancePairCAdd. m_result is
// flattened n_dep x n_arr x n_add, column-major (R array convention).
struct phastCAdd : public Worker{

  CGraph *m_gr;
  IVec m_dep;
  IVec m_arr;
  int n_add;
  DVec m_result;

  phastCAdd(CGraph *gr, IVec dep, IVec arr);

  void operator()(std::size_t begin, std::size_t end);

};

#endif
