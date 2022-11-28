#include <queue>
#include <string>
#include <stack>
#include <chrono>
#include "graph.h"
#include "bush.h"
#include <RcppParallel.h>
#include <Rcpp.h>

using namespace std;



Bush::Bush(){

}

Bush::Bush(Graph* ptr, int origin, IVec m_arr, DVec m_dem, Bush_vectors* bv_ptr, double num_tol){
  is_equilibrated = false;
  changed = false;
  arr= m_arr;
  dem = m_dem;
  gptr = ptr;
  bv = bv_ptr;
  root = origin;
  tol = num_tol;


  flow.resize(gptr->nbedge, 0.0);
  edges.resize(gptr->nbedge, 0);
  order.resize(gptr->nbnode, 0);


  DVec distances(gptr->nbnode, numeric_limits<double>::max());
  IVec parents(gptr->nbnode, -1);
  IVec links(gptr->nbnode, -1);
  distances[origin] = 0.0;

  PQ Q;
  Q.push(make_pair(origin, 0.0));

  while (!Q.empty()) {
    int v = Q.top().first;
    double w = Q.top().second;
    Q.pop();

    if (w <= distances[v]) {

      for (unsigned int i = gptr->indG[v]; i< gptr->indG[v + 1]; i++) {

        int v2 = gptr->nodeG[i];
        double w2 = gptr->wG[i];

        if (distances[v] + w2 < distances[v2]) {
          distances[v2] = distances[v] + w2;
          parents[v2] = v;
          links[v2] = i;
          Q.push(make_pair(v2, distances[v2]));
        }
      }
    }

  }

  for (unsigned int i = 0; i < gptr->nbnode; i++){
    if (links[i] == -1) continue;
    edges[links[i]] = 1;
  }

}

void Bush::ordering(){

  queue <int> sel;

  sel.push(root);

  // all nodes with 0 incoming edge
  fill(bv->incoming.begin(), bv->incoming.end(), 0);
  for (int i = 0; i < gptr->nodeG.size(); i++){
    if (edges[i] == 1) bv->incoming[gptr->nodeG[i]] += 1;
  }


  // ordering

  int count = 0;
  while (!sel.empty()){
    int v = sel.front();


    sel.pop();
    order[count] = v;
    count += 1;
    for (int i = gptr->indG[v]; i < gptr->indG[v + 1]; i++){

      if (edges[i] == 1){
        bv->incoming[gptr->nodeG[i]] += -1;
        if (bv->incoming[gptr->nodeG[i]] == 0) {
          sel.push(gptr->nodeG[i]);

        }
      }

    }


  }

}


void Bush::minmaxtree2(int mode){

  bv->reinitialize();
  bv->sdist[root] = 0.0;
  bv->ldist[root] = 0.0;

  for (int i = 0; i < edges.size(); i++){
    if (edges[i] == 1){
      bv->incoming[gptr->nodeG[i]] += 1;
    }

  }

  int count = 0;

  while (count < gptr->nbnode){
    int v = order[count];

    // shortest path
    for (int i = gptr->indGr[v]; i < gptr->indGr[v + 1]; i++){

      int w = gptr->nodeGr[i];
      int is_edge = 0;
      int edge_index = 0;


      for (int j = gptr->indG[w]; j < gptr->indG[w + 1]; j++){

        double ww = gptr->wG[j];

        if (gptr->nodeG[j] == v && edges[j] == 1){
          is_edge = edges[j];
          edge_index = j;
          if (is_edge == 1 && bv->sdist[v] > bv->sdist[w] + ww) {
            bv->sdist[v] = bv->sdist[w] + ww;
            bv->sparents[v] = w;
            bv->slink[v] = j;

          }

        }
      }

    }

    //longest path
    for (int i = gptr->indGr[v]; i < gptr->indGr[v + 1]; i++){

      int w = gptr->nodeGr[i];
      int is_edge = 0;
      int edge_index = 0;


      for (int j = gptr->indG[w]; j < gptr->indG[w + 1]; j++){

        double ww = gptr->wG[j];

        if (gptr->nodeG[j] == v && edges[j] == 1){
          is_edge = edges[j];
          edge_index = j;
          if (bv->incoming[v] > 1 && is_edge == 1 && bv->ldist[v] < bv->ldist[w] + ww && (mode== 0 || (mode == 1 && (flow[edge_index] > tol || (bv->slink[v] == j))) || (mode==2 && flow[edge_index] > tol))) {

            bv->ldist[v] = bv->ldist[w] + ww;
            bv->lparents[v] = w;
            bv->llink[v] = j;
          }

          if (bv->incoming[v] == 1){
            bv->ldist[v] = bv->ldist[w] + ww;
            bv->lparents[v] = w;
            bv->llink[v] = j;
          }
        }
      }

    }


    count += 1;
  }


}

void Bush::loadAON(){

  for (int i = 0; i < arr.size(); i++){
    int node = arr[i];

    // the other ones
    for (auto p = node; p != -1; p = bv->sparents[p]){
      int prev = bv->sparents[p];

      if (prev == -1) break;

      for (int j = gptr->indG[prev]; j < gptr->indG[prev + 1]; j++){
        if (gptr->nodeG[j] == p && edges[j] == 1 && bv->slink[p] == j){
          flow[j] += dem[i];
          gptr->flow[j] += dem[i];
          break;
        }
      }
    }

  }
}


void Bush::optimize2(){

  changed = false;
  // add shortcuts
  IVec added(edges.size(), 0);

  for (int i = 0; i < edges.size(); i++){
    int w = gptr->indG2[i];
    double ww = gptr->wG[i];
    int v = gptr->nodeG[i];


    int was_in_bush = edges[i];
    if (flow[i] < tol) flow[i] = 0.0;

    // edges[i] = 0;

    if (bv->sdist[v] == numeric_limits<double>::max()) continue;
    if (flow[i] > 0.0) {
      edges[i] = 1;
      continue;
    }

    if (bv->ldist[w] == -numeric_limits<double>::max() && bv->ldist[v] > -numeric_limits<double>::max()) continue;

    if (was_in_bush == 0 && (bv->ldist[v] > bv->ldist[w]) && (bv->sdist[v] > bv->sdist[w] + ww)) {
      //  Rcpp::Rcout << "crit1 " << w << " -> " << v  << " sp1:" <<sdist[w] << " sp2: " << sdist[v] << "lp1: " << ldist[w] << " lp2: " << ldist[v] << endl;
      edges[i] = 1;
      added[i] = 1;
      changed = true;
    }else if (was_in_bush == 1 && flow[i] == 0.0 &&  bv->sparents[v] == w && bv->slink[v] == i){ //   &&  sdist[v] == sdist[w] + ww
      //  Rcpp::Rcout << "crit2 " << w << " -> " << v  << " sp1:" <<sdist[w] << " sp2: " << sdist[v] << "lp1: " << ldist[w] << " lp2: " << ldist[v] << "edge: " << gptr->wG[i] << endl;
      edges[i] = 1;
      added[i] = 1;
    }
  }

  //Rcpp::Rcout << root << " " << newarcs << endl;


  fill(bv->incoming.begin(), bv->incoming.end(), 0);
  for (int i = 0; i < edges.size(); i++){
    if (edges[i] == 1){
      bv->incoming[gptr->nodeG[i]] += 1;
    }

  }

  // remove zero flow links

  for (int i = 0; i < edges.size(); i++){

    if (edges[i] == 1 && flow[i] == 0.0 && bv->incoming[gptr->nodeG[i]] > 1 &&  added[i] == 0) {// && added[i] == 0
      edges[i] = 0;
      bv->incoming[gptr->nodeG[i]] -= 1;
      changed = true;
      // Rcpp::Rcout<< "rm : " << gptr->indG2[i] << " -> " << gptr->nodeG[i] << endl;
    }

  }



}



void Bush::update_one_cost(int edge_index){
  double old_w = gptr->wG[edge_index];

  gptr->wG[edge_index] = gptr->ftt[edge_index] * (1.0 + gptr->alpha[edge_index] * (pow(gptr->flow[edge_index]/gptr->cap[edge_index], gptr->beta[edge_index])));
  // reversed
  for (int j = gptr->indGr[gptr->nodeG[edge_index]]; j < gptr->indGr[gptr->nodeG[edge_index] + 1]; j++){
    if (gptr->nodeGr[j] == gptr->indG2[edge_index] && gptr->wGr[j] == old_w){
      gptr->wGr[j] = gptr->wG[edge_index];
      break;
    }
  }

  // update data structure
  for (int j = 0; j < gptr->data[gptr->indG2[edge_index]].size(); j++){
    if (gptr->data[gptr->indG2[edge_index]][j].first == gptr->nodeG[edge_index] && gptr->data[gptr->indG2[edge_index]][j].second == old_w){
      gptr->data[gptr->indG2[edge_index]][j].second = gptr->wG[edge_index];
      break;
    }
  }
}

void Bush::shift_flow(int node){

  // find pair of alternative segments (PAS)
  sedges.resize(0);
  ledges.resize(0);
  sdiff.resize(0);
  ldiff.resize(0);

  // divergence node
  // first, set of shortest path node
  int div = -1;
  set<int> sp_nodes;

  //sp_nodes.insert(node);
  for (auto p = node; p != -1; p = bv->sparents[p]) {
    if (bv->sparents[p] == -1) break;
    int from = bv->sparents[p];
    sp_nodes.insert(from);

  }

  for (auto p = node; p != -1; p = bv->lparents[p]) {
    if (bv->lparents[p] == -1) break;
    int from = bv->lparents[p];



    if (sp_nodes.find(from) != sp_nodes.end()){

      div = from;
      break;
    }
  }



  // shortest path edge indexes
  for (auto p = node; p != -1; p = bv->sparents[p]) {
    if (bv->sparents[p] == -1) continue;
    int from = bv->sparents[p];


    for (int i = gptr->indG[from]; i < gptr->indG[from + 1]; i++){
      if (gptr->nodeG[i] == p && edges[i] == 1 && bv->slink[p] == i){
        sedges.push_back(i);
        break;
      }
    }

    if (from == div) break;

  }

  if (div == -1) return;

  // longest path edge indexes
  for (auto p = node; p != -1; p = bv->lparents[p]) {
    if (bv->lparents[p] == -1) continue;
    int from = bv->lparents[p];



    for (int i = gptr->indG[from]; i < gptr->indG[from + 1]; i++){
      if (gptr->nodeG[i] == p && edges[i] == 1 && bv->llink[p] == i){
        ledges.push_back(i);
        break;
      }
    }

    if (from == div) break;

  }


  // identify divergent edges

  sort(sedges.begin(), sedges.end());
  sort(ledges.begin(), ledges.end());
  set_difference(sedges.begin(), sedges.end(), ledges.begin(), ledges.end(), std::back_inserter(sdiff));
  set_difference(ledges.begin(), ledges.end(), sedges.begin(), sedges.end(), std::back_inserter(ldiff));

  // newton shift
  // cost difference
  double num = 0.0;

  for (int i = 0; i < ldiff.size(); i++){
    // if (root == 12 && node == 14) Rcpp::Rcout << "long "<< gptr->indG2[ldiff[i]] << " -> " << gptr->nodeG[ldiff[i]] << " : " << gptr->wG[ldiff[i]] << endl;
    num += gptr->wG[ldiff[i]];
  }

  for (int i = 0; i < sdiff.size(); i++){
    //  if (root == 12 && node == 14) Rcpp::Rcout << "short "<< gptr->indG2[sdiff[i]] << " -> " << gptr->nodeG[sdiff[i]] << " : " << gptr->wG[sdiff[i]] << endl;
    num -= gptr->wG[sdiff[i]];
  }



  // sum of derivatives
  double den = 0.0;
  double minflow2 = numeric_limits<double>::max();
  for (int j = 0; j < sdiff.size(); j++){
    int i = sdiff[j];

    if (flow[i] < minflow2) minflow2 = flow[i];

    den += gptr->ftt[i] * (gptr->alpha[i] * gptr->beta[i] * pow(gptr->flow[i]/gptr->cap[i], gptr->beta[i] - 1) )/gptr->cap[i];

  }

  double minflow = numeric_limits<double>::max();
  for (int j = 0; j < ldiff.size(); j++){
    int i = ldiff[j];

    if (flow[i] < minflow) minflow = flow[i];

    den += gptr->ftt[i] * (gptr->alpha[i] * gptr->beta[i] * pow(gptr->flow[i]/gptr->cap[i], gptr->beta[i] - 1) )/gptr->cap[i];

  }

  if (den == 0.0) den = 0.000001;
  double shift = min(minflow, num/den);
  shift = max(shift, -minflow2);
  //if (shift< 0) shift=0;



  // equilibrate flows
  //  if (sdiff.size()> 0)  Rcpp::Rcout<< "origin : " << root << " node:  " << node << " div: " << div << " num:" << num << " den:" << den << " shift:" << shift << " flow1:" << minflow << " flow2:" << minflow2 << endl;
  //  if (root == 0) Rcpp::Rcout << "flow 1er edge : " << flow[0] << endl;
  for (int j = 0; j < sdiff.size(); j++) {
    //if (abs(shift) < tol) continue;
    if (sdiff.size() > 0 && ldiff.size() >0){ //!same_last_edge && incoming[node] > 1

      // shift graph flow
      gptr->flow[sdiff[j]] += shift;
      // shift bush flow
      flow[sdiff[j]] += shift;

      if (gptr->flow[sdiff[j]] < 0.0) gptr->flow[sdiff[j]] = 0.0;
      if (flow[sdiff[j]] < 0) flow[sdiff[j]] = 0.0;

      //update cost
      update_one_cost(sdiff[j]);

    }

  }
  for (int j = 0; j < ldiff.size(); j++) {
    //if (abs(shift) < tol) continue;

    if (sdiff.size() > 0 && ldiff.size() >0) {
      gptr->flow[ldiff[j]] -= shift;
      flow[ldiff[j]] -= shift;

      if (gptr->flow[ldiff[j]] < 0.0) gptr->flow[ldiff[j]] = 0.0;
      if (flow[ldiff[j]] < 0) flow[ldiff[j]] = 0.0;

      // update cost
      update_one_cost(ldiff[j]);
    }
  }

}


