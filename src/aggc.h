#ifndef AGGC_H
#define AGGC_H

#include <vector>
#include <RcppParallel.h>
#include <Rcpp.h>
#include "graph.h"
#include "cgraph.h"
using namespace RcppParallel;

// aggregate additionnal weight from a normal graph to contracted graph
struct aggC : public Worker{
  CGraph* m_gr;
  Graph* m_or;
  DVec m_result; // upgraph
  DVec m_result2; // downgraph

  aggC(CGraph* gr, Graph* original);

  // iterate over m_gr->nbnode
  void operator()(std::size_t begin, std::size_t end);
};

// Same as aggC, generalized to N aux columns at once: for every shortcut
// edge in the contracted graph, sums each of the original graph's n_add
// columns over the unpacked original edges (instead of a single running
// sum), so distancePairCAdd never needs to unpack shortcuts during queries.
struct aggCAdd : public Worker{
  CGraph* m_gr;
  Graph* m_or;
  int n_add;
  DVec m_result; // flattened (m_gr->nodeG.size() x n_add), column-major, upgraph
  DVec m_result2; // flattened (m_gr->nodeGr.size() x n_add), column-major, downgraph

  aggCAdd(CGraph* gr, Graph* original);

  // iterate over m_gr->nbnode
  void operator()(std::size_t begin, std::size_t end);
};

#endif
