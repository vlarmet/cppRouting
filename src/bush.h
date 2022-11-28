#ifndef BUSH_H
#define BUSH_H

#include <queue>
#include <string>
#include "graph.h"
#include <RcppParallel.h>
#include <Rcpp.h>
using namespace std;

struct Bush_vectors {
  DVec sdist;
  DVec ldist;
  IVec sparents;
  IVec lparents;
  IVec slink;
  IVec llink;
  IVec incoming;

  Bush_vectors(int n_node){
    sdist.resize(n_node, numeric_limits<double>::max());
    ldist.resize(n_node, -numeric_limits<double>::max());
    sparents.resize(n_node, -1);
    lparents.resize(n_node, -1);
    slink.resize(n_node, -1);
    llink.resize(n_node, -1);
    incoming.resize(n_node, 0);
  }

  void reinitialize(){
    fill(sdist.begin(), sdist.end(), numeric_limits<double>::max());
    fill(ldist.begin(), ldist.end(), -numeric_limits<double>::max());
    fill(sparents.begin(),sparents.end(), -1);
    fill(lparents.begin(), lparents.end(), -1);
    fill(slink.begin(), slink.end(), -1);
    fill(llink.begin(), llink.end(), -1);
    fill(incoming.begin(), incoming.end(), 0);
  }
};


class Bush {
public:
  int root;
  bool is_equilibrated;
  bool changed;
  double tol;

  // for flow shifting
  IVec sedges;
  IVec ledges;
  IVec sdiff;
  IVec ldiff;

  // persistent vectors (stored for each bush)
  IVec edges;
  IVec order;
  DVec flow;
  /*
  DVec sdist;
  DVec ldist;
  IVec sparents;
  IVec lparents;
  IVec slink;
  IVec llink;
  IVec incoming;
   */
  // vectors shared between bushes, no need to store them individually
  Bush_vectors* bv;

  // graph object (address)
  Graph* gptr;

  // OD matrix
  IVec arr;
  DVec dem;

  Bush(Graph* ptr, int origin, IVec m_arr, DVec m_dem, Bush_vectors* bv_ptr, double num_tol);
  Bush();


  void minmaxtree2(int mode);
  void loadAON();

  void optimize2();
  void ordering();
  void shift_flow(int node);

  void update_one_cost(int edge_index);


};

#endif
