#ifndef PATH_PAIR_H
#define PATH_PAIR_H

#include <vector>
#include <RcppParallel.h>
#include <Rcpp.h>
#include "graph.h"
#include "cgraph.h"
using namespace RcppParallel;


struct pathPair : public Worker{

  // graph pointer to avoid copying input
  Graph *m_gr;
  IVec m_dep;
  IVec m_arr;
  IVec m_keep;
  double m_lim; // used by isochrone OR detour
  int algorithm; // 0 : dijkstra with early stop, 1 : bidirectional dijkstra, 2 : A*, 3 : NBA*
  vector<SVec> m_result;

  // constructor
  pathPair(Graph *gr, IVec dep, IVec arr, IVec keep, double lim, int algo);

  // algorithms as member functions

  void dijkstra_early_stop(std::size_t begin, std::size_t end);
  void bidir(std::size_t begin, std::size_t end);
  void astar(std::size_t begin, std::size_t end);
  void nba(std::size_t begin, std::size_t end);
  void iso(std::size_t begin, std::size_t end);
  void detour(std::size_t begin, std::size_t end);

  // overload operator ()
  void operator()(std::size_t begin, std::size_t end);

};


struct pathPairC : public Worker{

  // graph pointer to avoid copying input
  CGraph *m_gr;
  IVec m_dep;
  IVec m_arr;
  IVec m_keep;
  int algorithm;
  vector<SVec> m_result;

  // constructor
  pathPairC(CGraph *gr, IVec dep, IVec arr, IVec keep, int algo);

  // algorithms as member functions

  void bidirmod(std::size_t begin, std::size_t end);


  // overload operator ()
  void operator()(std::size_t begin, std::size_t end);

};


#endif
