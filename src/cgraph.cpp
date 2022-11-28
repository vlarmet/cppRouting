#include "graph.h"
#include "cgraph.h"
#include "distance_pair.h"
#include "path_pair.h"
#include "distance_mat.h"
#include "path_mat.h"
#include <vector>
#include <RcppParallel.h>
#include <Rcpp.h>
using namespace std;


// Constructor
/// from R, not contracted
CGraph::CGraph(IVec &gfrom, IVec &gto, DVec &gw, int nb)
{
  nbnode = nb;
  nbedge = gfrom.size();
  cdata.resize(nb);
  is_contracted = false;
  rank.resize(nb);
  add.resize(0);
  addr.resize(0);

  for (int i = 0; i != nbedge; ++i) {
    bool exist = false;
    for (int j = 0; j < cdata[gfrom[i]].size(); j++){
      if (cdata[gfrom[i]][j].first == gto[i]){
        exist = true;
        if (cdata[gfrom[i]][j].second > gw[i]){
          cdata[gfrom[i]][j].second = gw[i];
        }
        break;
      }
    }
    if (!exist) cdata[gfrom[i]].push_back(make_pair(gto[i], gw[i]));
  }

  fill(rank.begin(), rank.end(), 0);

}

/// from R, already contracted
CGraph::CGraph(IVec &gfrom, IVec &gto, DVec &gw, int nb, IVec &m_rank, IVec &m_shortf, IVec &m_shortt, IVec &m_shortc, bool phast)
{
  nbnode = nb;
  nbedge = gfrom.size();
  cdata.resize(nb);
  cdataR.resize(nb);
  is_contracted = true;
  rank = m_rank;
  shortf = m_shortf;
  shortt = m_shortt;
  shortc = m_shortc;

  add.resize(0);
  addr.resize(0);

  if (phast){
    for (int i = 0; i != nbedge; ++i) {
      if(gfrom[i] > gto[i]) cdata[gfrom[i]].push_back(make_pair(gto[i], gw[i]));
      if(gfrom[i] < gto[i]) cdataR[gto[i]].push_back(make_pair(gfrom[i], gw[i]));
    }

  } else{
    for (int i = 0; i != nbedge; ++i) {
      if(rank[gfrom[i]] < rank[gto[i]]) cdata[gfrom[i]].push_back(make_pair(gto[i], gw[i]));
      if(rank[gfrom[i]] > rank[gto[i]]) cdataR[gto[i]].push_back(make_pair(gfrom[i], gw[i]));
    }

  }


}

// From a graph object (c++), not contracted
CGraph::CGraph(Graph* graph):
  nbnode(graph->nbnode), nbedge(graph->nbedge),cdata(graph->nbnode), is_contracted(false), rank(graph->nbnode)
{
  add.resize(0);
  addr.resize(0);
  //cdata = graph->data;
  cdata.resize(graph->data.size());

  int count= 0;
  for (int i = 0; i != graph->data.size(); ++i) {

    for (int k = 0; k < graph->data[i].size(); k++){
      bool exist = false;

      for (int j = 0; j < cdata[i].size(); j++){
        if (cdata[i][j].first == graph->data[i][k].first){
          exist = true;
          if (cdata[i][j].second > graph->data[i][k].second){
            cdata[i][j].second = graph->data[i][k].second;
          }
          break;
        }
      }
      if (!exist) {
        cdata[i].push_back(graph->data[i][k]);
        count ++ ;
      }

    }
  }

  nbedge = count;

  fill(rank.begin(), rank.end(), 0);
}


G CGraph::getReverse(G &graph){
  G dataR(nbnode);
  for (int i = 0; i != nbnode; ++i) {

    for (int j = 0; j < graph[i].size(); ++j){
      int v = graph[i][j].first;
      double w = graph[i][j].second;

      dataR[v].push_back(make_pair(i, w));

    }

  }
  return dataR;
}

void CGraph::setDict(SVec &m_dict){
  dict = m_dict;
}

void CGraph::construct_shortcuts(){

  shortcuts.resize(nbnode);

  for (int i=0; i < shortf.size(); i++){

    shortcuts[shortf[i]].push_back(make_pair(shortt[i], shortc[i]));

  }
}

Rcpp::List CGraph::getEdges(){
  Rcpp::List edges(3);

  IVec Newfrom(nbedge);
  IVec Newto(nbedge);
  DVec Neww(nbedge);

  int index=0;
  for (int i=0; i < nbnode;i++){
    for (int j=0; j < cdata[i].size();j++){

      Newfrom[index]=i;
      Newto[index]= cdata[i][j].first;
      Neww[index]= cdata[i][j].second;
      index+=1;

    }
  }

  edges[0] = Newfrom;
  edges[1] = Newto;
  edges[2] = Neww;

  return edges;
}

/*
void CGraph::unpack(IVec &path){
  IVec temp(path);

  while (true){


    int again=0;
    vector<pair<int,int> > to_insert;
    int count=0;
    for (int i=0; i < (path.size() -1);i++){

      for (int j = 0; j < shortcuts[path[i]].size(); j++){

        if (shortcuts[path[i]][j].first == path[i+1]){
          to_insert.push_back(make_pair(i+count+1, shortcuts[path[i]][j].second));
          count+=1;
          again=1;
          break;
        }
      }

    }
    for (int i=0; i < count; i++) temp.insert(temp.begin() + to_insert[i].first, to_insert[i].second);

    path = temp;

    if (again==0) break;

  }


}
 */

void CGraph::unpack(IVec &path){

  while (true){


    int again=0;
    vector<pair<int,int> > to_insert;
    int count=0;
    for (int i=0; i < (path.size() -1);i++){

      for (int j = 0; j < shortcuts[path[i]].size(); j++){

        if (shortcuts[path[i]][j].first == path[i+1]){
          to_insert.push_back(make_pair(i+count+1, shortcuts[path[i]][j].second));
          count+=1;
          again = 1;

          break;
        }

      }


    }
    for (int i=0; i < count; i++) path.insert(path.begin() + to_insert[i].first, to_insert[i].second);

    if (again==0) break;

  }


}

void CGraph::to_adj_list(bool reversed, bool sorting){
  int count = 0;
  if (!reversed){
    for (int i = 0; i < nbnode; i++){
      count += cdata[i].size();
    }

    nodeG.resize(count);
    wG.resize(count);
    indG.resize(nbnode + 1);
    count=0;

    if (sorting){
      for (int i = 0; i < nbnode; i++) sort(cdata[i].begin(), cdata[i].end());
    }

    for (int i=0; i < cdata.size();i++){
      indG[i] = count;

      for (int j = 0; j < cdata[i].size(); j++){
        nodeG[count] = cdata[i][j].first;
        wG[count] = cdata[i][j].second;
        count += 1;
      }
    }
    indG[nbnode]=count;

  } else{

    for (int i = 0; i < nbnode; i++){
      count += cdataR[i].size();
    }

    nodeGr.resize(count);
    wGr.resize(count);
    indGr.resize(nbnode + 1);
    count=0;

    if (sorting){
      for (int i = 0; i < nbnode; i++) sort(cdataR[i].begin(), cdataR[i].end());
    }

    for (int i=0; i < cdataR.size();i++){
      indGr[i] = count;

      for (int j = 0; j < cdataR[i].size(); j++){
        nodeGr[count] = cdataR[i][j].first;
        wGr[count] = cdataR[i][j].second;
        count += 1;
      }
    }
    indGr[nbnode] = count;
  }

}


///////////////////////////////////////////// Algorithms
bool CGraph::stall(int &node, DVec &distances, G &graph){
  bool res = false;
  for (int i = 0; i < graph[node].size(); i++){

    if ((distances[graph[node][i].first] + graph[node][i].second) < distances[node]){

      res = true;
      break;
    }

  }

  return res;
}


DVec CGraph::routing_dvec(IVec dep, IVec arr, int algo){
  distancePairC dijfunc(this, dep, arr, algo);
  parallelFor(0, dep.size(), dijfunc, 1, 12);
  return dijfunc.m_result;
}

vector<SVec> CGraph::routing_svec(IVec dep, IVec arr, IVec keep, int algo){
  pathPairC dijfunc(this, dep, arr, keep, algo);
  parallelFor(0, dep.size(), dijfunc, 1, 12);
  return dijfunc.m_result;
}

Rcpp::NumericMatrix CGraph::routing_dmat(IVec dep, IVec arr, int algo){
  Rcpp::NumericMatrix final_result;
  // buckets
  if (algo == 0){

    G BucketF(nbnode);

    //Forward

    DVec distances(nbnode, numeric_limits<double>::max());
    IVec visited(nbnode, 0);


    for (unsigned int k=0; k!=dep.size();k++){
      int start =dep[k];
      distances[start] = 0.0;
      PQ Q;
      Q.push(make_pair(start, 0.0));

      while (true) {

        if (Q.empty()){
          break;
        }

        if (!Q.empty()){

          int v=Q.top().first;
          double w=Q.top().second;
          Q.pop();
          visited[v] = 1;


          if (w <= distances[v] && !stall(v, distances, cdataR)) {
            for (int i=0; i< cdata[v].size(); i++){
              int v2 = cdata[v][i].first;
              double w2 = cdata[v][i].second;


              if (distances[v] + w2 < distances[v2]) {
                distances[v2] = distances[v] + w2;
                Q.push(make_pair(v2, distances[v2]));
                visited[v2] = 1;

              }
            }
          }

        }


      }

      for (int i=0; i < distances.size();i++){
        if (distances[i] < numeric_limits<double>::max()){
          BucketF[i].push_back(make_pair(k, distances[i]));
        }

      }


      fill(distances.begin(), distances.end(), numeric_limits<double>::max());


    }

    int count2=0;
    for (int i=0; i < BucketF.size();i++){
      count2+= BucketF[i].size();
      std::sort(BucketF[i].begin(),BucketF[i].end());
    }


    //Backward in parallel
    //Graph vectors
    to_adj_list(true, false);
    cdata = BucketF;
    to_adj_list(false, false);

    G ().swap(BucketF);

    Rcpp::NumericMatrix result(arr.size(),dep.size());


    distanceMatC Dijfunc(this, dep, arr, result);
    parallelFor(0,arr.size(),Dijfunc);


    final_result = Rcpp::transpose(result);
  }

  // phast
  if (algo == 1){

    Rcpp::NumericMatrix result(dep.size(),arr.size());

    to_adj_list(false, true);
    to_adj_list(true, true);

    phastC Dijfunc(this, dep, arr,result);

    parallelFor(0,dep.size(),Dijfunc);
    final_result = result;
  }

  return final_result;
}

vector<SVec> CGraph::routing_smat(IVec dep, IVec arr, IVec keep){
  pathMatC dijfunc(this, dep, arr, keep);
  parallelFor(0, dep.size(), dijfunc);
  return dijfunc.m_result;
}

