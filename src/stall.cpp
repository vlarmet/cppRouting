#include "stall.h"
#include "graph.h"

bool Stall_par(int &node,
               DVec &distance,
               IVec &NodeG,
               DVec &WG,
               IVec &IndG){

  bool res = false;
  for (int i = IndG[node]; i < IndG[node + 1]; i++){
    if ((distance[NodeG[i]] + WG[i]) < distance[node]){
      res = true;
      break;
    }
  }

  return res;
}
