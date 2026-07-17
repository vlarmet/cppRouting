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

// Distance-only pair query that sums N aux columns inline during the same
// search (no path reconstruction), generalizing distancePair's single-column
// add/addr mechanism to n_add columns at once.
struct distancePairAdd : public Worker{

  Graph *m_gr;
  IVec m_dep;
  IVec m_arr;
  int n_add;
  int algorithm; // 0 : dijkstra with early stop, 1 : bidirectional dijkstra, 2 : A*, 3 : NBA*
  DVec m_result; // flattened n_pairs x n_add, column-major

  distancePairAdd(Graph *gr, IVec dep, IVec arr, int algo);

  void dijkstra_early_stop(std::size_t begin, std::size_t end);
  void bidir(std::size_t begin, std::size_t end);
  void astar(std::size_t begin, std::size_t end);
  void nba(std::size_t begin, std::size_t end);

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

// Same as distancePairC (contracted graph, modified bidirectional dijkstra),
// generalized to n_add aux columns summed inline. Aux values on m_gr->add /
// m_gr->addr are per-shortcut sums already precomputed by aggCAdd, so no
// unpacking happens during the query itself.
struct distancePairCAdd : public Worker{

  CGraph *m_gr;
  IVec m_dep;
  IVec m_arr;
  int n_add;
  DVec m_result; // flattened n_pairs x n_add, column-major

  distancePairCAdd(CGraph *gr, IVec dep, IVec arr);

  void bidirmod(std::size_t begin, std::size_t end);

  void operator()(std::size_t begin, std::size_t end);

};


#endif
