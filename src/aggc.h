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

#endif
