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

#endif
