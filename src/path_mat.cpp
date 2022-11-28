// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]

#include <vector>
#include <Rcpp.h>
#include <RcppParallel.h>
#include "graph.h"
#include "path_mat.h"
using namespace RcppParallel;

// constructor
pathMat::pathMat(Graph *gr, IVec dep, IVec arr, IVec keep, DVec lim, bool setdif, int algo) :
  m_gr(gr), m_dep(dep), m_arr(arr), m_keep(keep), m_lim(lim), m_setdif(setdif), algorithm(algo)
{
  m_result.resize(m_dep.size());

}

// switch between algorigthms
void pathMat::operator()(std::size_t begin, std::size_t end){
  if (algorithm == 0) dijkstra_path_mat(begin, end);
  if (algorithm == 1) multi_iso(begin, end);

}


// dijkstra multi-paths

void pathMat::dijkstra_path_mat(std::size_t begin, std::size_t end){
  DVec distances(m_gr->nbnode, numeric_limits<double>::max());
  IVec parents(m_gr->nbnode, -1);

  int count = 0;

  for (std::size_t k = begin; k != end; k++){

    int start = m_dep[k];

    distances[start] = 0.0;

    PQ Q;
    Q.push(make_pair(start, 0.0));

    while (!Q.empty()) {
      int v = Q.top().first;
      double w = Q.top().second;
      Q.pop();

      if (w <= distances[v]) {

        for (int i = m_gr->indG[v]; i<  m_gr->indG[v+1]; i++) {

          int v2 = m_gr->nodeG[i];
          double w2 = m_gr->wG[i];

          if (distances[v] + w2 < distances[v2]) {
            distances[v2] = distances[v] + w2;
            parents[v2] = v;
            Q.push(make_pair(v2, distances[v2]));
          }
        }
      }

    }



    SVec result(m_arr.size());


    for (unsigned int i=0; i < m_arr.size();i++){
      int end = m_arr[i];
      SVec result2;
      string s;


      for (auto p = end; p != -1; p = parents[p]){
        if (m_keep[p] == 1) {

          s.append(m_gr->dict[p]);
          s.append(" ");
        }
      }

      result[i] = s;

    }


    m_result[k] = result;

    //Reinitialize
    fill(distances.begin(),distances.end(),numeric_limits<double>::max());
    fill(parents.begin(),parents.end(),-1);


  }


}


// multi isochrones

void pathMat::multi_iso(std::size_t begin, std::size_t end){

  double max_limit = *max_element(m_lim.begin(), m_lim.end());
  DVec distances(m_gr->nbnode, numeric_limits<double>::max());

  for (std::size_t k = begin; k != end; k++){

    int start = m_dep[k];
    distances[start] = 0.0;
    PQ Q;
    Q.push(make_pair(start, 0.0));

    while (!Q.empty()) {
      int v = Q.top().first;
      double w = Q.top().second;
      Q.pop();

      if (w <= distances[v]) {
        for (int i = m_gr->indG[v]; i <  m_gr->indG[v+1]; i++) {

          int v2 = m_gr->nodeG[i];
          double w2 = m_gr->wG[i];

          if (distances[v] + w2 < distances[v2]) {
            distances[v2] = distances[v] + w2;

            Q.push(make_pair(v2, distances[v2]));
          }

        }

      }
      if (distances[v] > max_limit){
        break;
      }
    }

    SVec result(m_lim.size());


    for (int i = 0; i < m_lim.size(); i++){
      string s;
      for (int j = 0; j < distances.size(); j++){
        if (m_keep[j] == 1){
          if (distances[j] < m_lim[i]){
            if (m_setdif){
              distances[j] = numeric_limits<double>::max();
            }

           // result[i].push_back(m_gr->dict[j]);
            s.append(m_gr->dict[j]);
            s.append(" ");
          }
        }
      }
      result[i] = s;
    }

    m_result[k] = result;
    //Reinitialize
    fill(distances.begin(),distances.end(),numeric_limits<double>::max());

  }

}

/////////////////////////////////////////////// contracted graph

pathMatC::pathMatC(CGraph *gr, IVec dep, IVec arr, IVec keep) :
  m_gr(gr), m_dep(dep), m_arr(arr), m_keep(keep)
{
  m_result.resize(m_dep.size());
}


void pathMatC::operator()(std::size_t begin, std::size_t end){

  DVec Distances(m_gr->nbnode, std::numeric_limits<double>::max());
  IVec parents(m_gr->nbnode, -1);
  IVec parents2(m_gr->nbnode, -1);

  for (std::size_t k=begin; k!=end;k++){

    int StartNode=m_dep[k];

    Distances[StartNode] = 0.0;
    PQ Q;
    Q.push(std::make_pair(StartNode, 0.0));


    while (true) {

      if (Q.empty()){
        break;
      }

      if (!Q.empty()){

        int v=Q.top().first;
        double w=Q.top().second;
        Q.pop();

        if (w <= Distances[v]) {
          for (int i = m_gr->indG[v]; i< m_gr->indG[v+1]; i++){
            int v2 = m_gr->nodeG[i];
            double w2 = m_gr->wG[i];


            if (Distances[v] + w2 < Distances[v2]) {
              Distances[v2] = Distances[v] + w2;
              parents[v2] = v;
              Q.push(std::make_pair(v2, Distances[v2]));

            }
          }
        }

      }
    }

    //Backward

    for (int i=0; i < (m_gr->indGr.size()-1); i++){

      for (int j=m_gr->indGr[i]; j < m_gr->indGr[i+1]; j++){


        if (Distances[m_gr->nodeGr[j]]+m_gr->wGr[j] < Distances[i]) {
          Distances[i] = Distances[m_gr->nodeGr[j]]+m_gr->wGr[j];
          parents2[i] = m_gr->nodeGr[j];
        }

      }
    }

    SVec result(m_arr.size());
    for (unsigned int i=0; i < m_arr.size();i++){
      int end = m_arr[i];

      IVec result2;
      int last_node = -1;
      for (auto p = end; p != -1; p = parents2[p]){
          result2.push_back(p);
          last_node = p;
      }


      reverse(result2.begin(),result2.end());

      if (last_node != -1){
        for (auto p = parents[last_node]; p != -1; p = parents[p]){
          result2.insert(result2.begin(), p);
          last_node = p;
        }
      }

      m_gr->unpack(result2);

      string s;
      for (int j = 0; j < result2.size(); j++){
        if (m_keep[result2[j]]){
          s.append(m_gr->dict[result2[j]]);
          s.append(" ");
        }

      }

      result[i] = s;

    }

    m_result[k] = result;

    std::fill(Distances.begin(),Distances.end(),std::numeric_limits<double>::max());
    std::fill(parents.begin(),parents.end(), -1);
    std::fill(parents2.begin(),parents2.end(), -1);

  }
}
