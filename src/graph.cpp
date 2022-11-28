#include "graph.h"
#include "distance_pair.h"
#include "path_pair.h"
#include "distance_mat.h"
#include "path_mat.h"
#include "aon.h"
#include <string>
#include <RcppParallel.h>
#include <Rcpp.h>
using namespace std;

// Constructor
Graph::Graph(vector<int> &gfrom, vector<int> &gto, vector<double> &gw, int nb){
  nbedge = gfrom.size();
  nbnode = nb;
  data = G(nbnode);
  add.resize(0);
  addr.resize(0);
  for (int i = 0; i != nbedge; ++i) {
    data[gfrom[i]].push_back(std::make_pair(gto[i], gw[i]));
  }

}

Graph::Graph(IVec &gfrom, IVec &gto, DVec &gw, int nb, bool rm_dup){
  nbedge = gfrom.size();
  nbnode = nb;
  data = G(nbnode);
  add.resize(0);
  addr.resize(0);

  for (int i = 0; i != nbedge; ++i) {
    bool exist = false;
    for (int j = 0; j < data[gfrom[i]].size(); j++){
      if (data[gfrom[i]][j].first == gto[i]){
        exist = true;
        if (data[gfrom[i]][j].second > gw[i]){
          data[gfrom[i]][j].second = gw[i];
        }
        break;
      }
    }
    if (!exist) data[gfrom[i]].push_back(make_pair(gto[i], gw[i]));
  }
}

Graph::Graph(IVec &gfrom, IVec &gto, DVec &gw, DVec &gadd, int nb){
  nbedge = gfrom.size();
  nbnode = nb;
  data = G(nbnode);
  vector<IVec> data2(nbnode);
  vector<IVec> data3(nbnode); // reversed
  for (int i = 0; i != nbedge; ++i) {
    data[gfrom[i]].push_back(std::make_pair(gto[i], gw[i]));
    data2[gfrom[i]].push_back(i);
    data3[gto[i]].push_back(i);
  }


  // adj list
  int count = 0;
  for (int i = 0; i < nbnode; i++){
    count += data[i].size();
  }

  nodeG.resize(count);
  wG.resize(count);
  add.resize(count);
  indG.resize(nbnode + 1);
  indG2.resize(count);

  count=0;

  for (int i=0; i < data2.size();i++){
    indG[i] = count;

    for (int j = 0; j < data2[i].size(); j++){
      nodeG[count] = gto[data2[i][j]];
      wG[count] = gw[data2[i][j]];
      add[count] = gadd[data2[i][j]];

      indG2[count] = i;
      count += 1;
    }
  }
  indG[nbnode]=count;

  // reversed

  nodeGr.resize(nodeG.size());
  wGr.resize(nodeG.size());
  addr.resize(nodeG.size());
  indGr.resize(nbnode + 1);

  count=0;

  for (int i=0; i < data3.size();i++){
    indGr[i] = count;

    for (int j = 0; j < data3[i].size(); j++){
      nodeGr[count] = gfrom[data3[i][j]];
      wGr[count] = gw[data3[i][j]];
      addr[count] = gadd[data3[i][j]];

      count += 1;
    }
  }
  indGr[nbnode]=count;

}

Graph::Graph(IVec &gfrom, IVec &gto, DVec &gw,
             DVec &gflow, DVec &gaux, DVec &gftt,
             DVec &galpha, DVec &gbeta, DVec &gcap, int nb){
  nbedge = gfrom.size();
  nbnode = nb;
  add.resize(0);
  addr.resize(0);

  data = G(nbnode);
  vector<IVec> data2(nbnode);
  for (int i = 0; i != nbedge; ++i) {
    data[gfrom[i]].push_back(std::make_pair(gto[i], gw[i]));
    data2[gfrom[i]].push_back(i);
  }

  // adj list
  int count = 0;
  for (int i = 0; i < nbnode; i++){
    count += data[i].size();
  }

  nodeG.resize(count);
  wG.resize(count);
  flow.resize(count);
  aux.resize(count);
  ftt.resize(count);
  alpha.resize(count);
  beta.resize(count);
  cap.resize(count);
  indG.resize(nbnode + 1);
  indG2.resize(count);

  count=0;

  for (int i=0; i < data2.size();i++){
    indG[i] = count;

    for (int j = 0; j < data2[i].size(); j++){
      nodeG[count] = gto[data2[i][j]];
      wG[count] = gw[data2[i][j]];
      flow[count] = gflow[data2[i][j]];
      aux[count] = gaux[data2[i][j]];
      ftt[count] = gftt[data2[i][j]];
      alpha[count] = galpha[data2[i][j]];
      beta[count] = gbeta[data2[i][j]];
      cap[count] = gcap[data2[i][j]];
      indG2[count] = i;
      count += 1;
    }
  }
  indG[nbnode]=count;


}

// Setters
void Graph::setReverse(){
  dataR = G(nbnode);
  for (int i = 0; i != nbnode; ++i) {

    for (int j = 0; j < data[i].size(); ++j){
      int v = data[i][j].first;
      double w = data[i][j].second;

      dataR[v].push_back(make_pair(i, w));

    }

  }
}

void Graph::setLatLon(DVec &m_lat, DVec &m_lon){
  lat = m_lat;
  lon = m_lon;
}

void Graph::setDict(SVec &m_dict){
  dict = m_dict;
}

Rcpp::List Graph::getEdges(){
  Rcpp::List edges(3);

  nbedge = 0;
  for (int i=0; i < data.size(); i++){
    nbedge += data[i].size();
  }

  IVec Newfrom(nbedge);
  IVec Newto(nbedge);
  DVec Neww(nbedge);

  int index=0;
  for (int i=0; i < data.size();i++){
    for (int j=0; j < data[i].size();j++){

      Newfrom[index]=i;
      Newto[index]= data[i][j].first;
      Neww[index]= data[i][j].second;
      index+=1;

    }
  }

  edges[0] = Newfrom;
  edges[1] = Newto;
  edges[2] = Neww;

  return edges;
}


void Graph::to_adj_list(bool reversed){
  int count = 0;
  if (!reversed){
    for (int i = 0; i < nbnode; i++){
      count += data[i].size();
    }

    nodeG.resize(count);
    wG.resize(count);
    indG.resize(nbnode + 1);
    count=0;

    for (int i=0; i < data.size();i++){
      indG[i] = count;

      for (int j = 0; j < data[i].size(); j++){
        nodeG[count] = data[i][j].first;
        wG[count] = data[i][j].second;
        count += 1;
      }
    }
    indG[nbnode]=count;

  } else{

    for (int i = 0; i < nbnode; i++){
      count += dataR[i].size();
    }

    nodeGr.resize(count);
    wGr.resize(count);
    if (add.size() > 0) addr.resize(count);
    indGr.resize(nbnode + 1);
    count=0;

    for (int i=0; i < dataR.size();i++){
      indGr[i] = count;

      for (int j = 0; j < dataR[i].size(); j++){
        nodeGr[count] = dataR[i][j].first;
        wGr[count] = dataR[i][j].second;
        count += 1;
      }
    }
    indGr[nbnode] = count;


  }
}


DVec Graph::routing_dvec(IVec dep, IVec arr, int algo){

  distancePair dijfunc(this, dep, arr, algo);
  parallelFor(0, dep.size(), dijfunc, 1, 12);
  return dijfunc.m_result;
}

vector<SVec> Graph::routing_svec(IVec dep, IVec arr, IVec keep, double lim, int algo){

  pathPair dijfunc(this, dep, arr, keep, lim, algo);
  parallelFor(0, dep.size(), dijfunc, 1, 12);
  return dijfunc.m_result;
}

Rcpp::NumericMatrix Graph::routing_dmat(IVec dep, IVec arr){
  Rcpp::NumericMatrix result(dep.size(), arr.size());
  distanceMat dijfunc(this, dep, arr, result);
  parallelFor(0, dep.size(), dijfunc, 1, 12);
  return result;
}

vector<SVec> Graph::routing_smat(IVec dep, IVec arr, IVec keep, DVec lim, bool setdif, int algo){

  pathMat dijfunc(this, dep, arr, keep, lim, setdif, algo);
  parallelFor(0, dep.size(), dijfunc, 1, 12);
  return dijfunc.m_result;
}
