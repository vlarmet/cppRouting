#ifndef PATH_MAT_H
#define PATH_MAT_H

#include <vector>
#include <RcppParallel.h>
#include <Rcpp.h>
#include "graph.h"
#include "cgraph.h"
using namespace RcppParallel;


struct pathMat : public Worker{

  // graph pointer to avoid copying input
  const Graph *m_gr;
  IVec m_dep;
  IVec m_arr;
  IVec m_keep;
  DVec m_lim;
  bool m_setdif;
  int algorithm; // 0 : multi path, 1 : isochrones with multiple limits

  //vector<vector<SVec>> m_result;
  vector<SVec> m_result;




  // constructor
  pathMat(Graph *gr, IVec dep, IVec arr, IVec keep, DVec lim, bool setdif, int algo);

  // algorithms as member functions

  void dijkstra_path_mat(std::size_t begin, std::size_t end);
  void multi_iso(std::size_t begin, std::size_t end);


  // overload operator ()
  void operator()(std::size_t begin, std::size_t end);

};


struct pathMatC : public Worker{

  // graph pointer to avoid copying input
  CGraph *m_gr;
  IVec m_dep;
  IVec m_arr;
  IVec m_keep;

  vector<SVec> m_result;




  // constructor
  pathMatC(CGraph *gr, IVec dep, IVec arr, IVec keep);


  // overload operator ()
  void operator()(std::size_t begin, std::size_t end);

};


#endif
