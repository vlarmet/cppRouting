#ifndef DISTANCE_PAIR_H
#define DISTANCE_PAIR_H

#include <vector>
#include <RcppParallel.h>
#include <Rcpp.h>
#include "graph.h"
#include "cgraph.h"
using namespace RcppParallel;


struct distancePair : public Worker{

  // graph pointer to avoid copying input
  Graph *m_gr;
  IVec m_dep;
  IVec m_arr;
  int algorithm; // 0 : dijkstra with early stop, 1 : bidirectional dijkstra, 2 : A*, 3 : NBA*
  bool add;
  DVec m_result;

  // constructor
  distancePair(Graph *gr, IVec dep, IVec arr, int algo);

  // algorithms as member functions

  void dijkstra_early_stop(std::size_t begin, std::size_t end);
  void bidir(std::size_t begin, std::size_t end);
  void astar(std::size_t begin, std::size_t end);
  void nba(std::size_t begin, std::size_t end);

  // overload operator ()
  void operator()(std::size_t begin, std::size_t end);

};

struct distancePairC : public Worker{

  // graph pointer to avoid copying input
  CGraph *m_gr;
  IVec m_dep;
  IVec m_arr;
  int algorithm;
  bool add;
  DVec m_result;

  // constructor
  distancePairC(CGraph *gr, IVec dep, IVec arr, int algo);

  // algorithms as member functions

  void bidirmod(std::size_t begin, std::size_t end);


  // overload operator ()
  void operator()(std::size_t begin, std::size_t end);

};


#endif
