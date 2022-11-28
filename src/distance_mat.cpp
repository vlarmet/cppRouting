// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]

#include <vector>
#include <Rcpp.h>
#include <RcppParallel.h>
#include "graph.h"
#include "cgraph.h"
#include "stall.h"
#include "distance_mat.h"
using namespace RcppParallel;

// constructor

distanceMat::distanceMat(Graph *gr, IVec dep, IVec arr, Rcpp::NumericMatrix result):
  m_gr(gr), m_dep(dep), m_arr(arr), m_result(result)
{
  add = false;
  if (m_gr->add.size() > 0) add = true;
}

void distanceMat::operator()(std::size_t begin, std::size_t end){
  dijkstra_mat(begin, end);

}

void distanceMat::dijkstra_mat(std::size_t begin, std::size_t end){

  DVec distances(m_gr->nbnode, numeric_limits<double>::max());
  DVec distances2;
  if (add) distances2.resize(m_gr->nbnode, numeric_limits<double>::max());

  for (std::size_t k = begin; k != end; k++){

    int start = m_dep[k];
    distances[start] = 0.0;
    if (add) distances2[start] = 0.0;

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
            if (add) distances2[v2] = distances2[v] + m_gr->add[i];

            Q.push(make_pair(v2, distances[v2]));
          }
        }
      }
    }


    for (int i = 0; i != m_arr.size(); i++){

      if (add){
        m_result(k,i)= distances2[m_arr[i]];
      } else{
        m_result(k,i)= distances[m_arr[i]];
      }



    }
    //Reinitialize
    fill(distances.begin(),distances.end(),numeric_limits<double>::max());
    if (add) fill(distances2.begin(),distances2.end(),numeric_limits<double>::max());
  }
}


///////////////////////////////////////////////////// contracted graph

distanceMatC::distanceMatC(CGraph *gr, IVec dep, IVec arr,
                           Rcpp::NumericMatrix result):
  m_gr(gr), m_dep(dep), m_arr(arr), m_result(result)
{

}


void distanceMatC::operator()(std::size_t begin, std::size_t end){

  DVec Distances(m_gr->nbnode, std::numeric_limits<double>::max());

  for (std::size_t k=begin; k!=end;k++){


    DVec Dist(m_dep.size(), std::numeric_limits<double>::max());

    int StartNode=m_arr[k];
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

        //Scan bucket

        if ((m_gr->indG[v+1] - m_gr->indG[v]) > 0){

          for (int i=m_gr->indG[v]; i < m_gr->indG[v+1];i++){

            int Source = m_gr->nodeG[i];
            double D = Distances[v] + m_gr->wG[i];

            if (Dist[Source] > D) {

              Dist[Source] = D;
            }

          }
        }

        if (w <= Distances[v]) {

          for (int i = m_gr->indGr[v]; i< m_gr->indGr[v+1]; i++){


            int v2 = m_gr->nodeGr[i];
            double w2 = m_gr->wGr[i];

            if (Distances[v] + w2 < Distances[v2]) {
              Distances[v2] = Distances[v] + w2;
              Q.push(std::make_pair(v2, Distances[v2]));

            }
          }
        }

      }


    }

    std::fill(Distances.begin(), Distances.end(), std::numeric_limits<double>::max());



    for (int i=0; i < Dist.size();i++){

      m_result(k,i)=Dist[i];

    }



  }
}


phastC::phastC(CGraph *gr, IVec dep, IVec arr,
               Rcpp::NumericMatrix result):
  m_gr(gr), m_dep(dep), m_arr(arr), m_result(result)
{
  add = false;
  if (m_gr->add.size() > 0) add = true;
}

void phastC::operator()(std::size_t begin, std::size_t end){

  DVec Distances(m_gr->nbnode, std::numeric_limits<double>::max());
  DVec distadd;
  if (add){
    distadd.resize(m_gr->nbnode, std::numeric_limits<double>::max());
  }

  for (std::size_t k=begin; k!=end;k++){

    int StartNode=m_dep[k];

    Distances[StartNode] = 0.0;
    if (add) distadd[StartNode] = 0.0;
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
          if (Stall_par(v, Distances, m_gr->nodeGr, m_gr->wGr, m_gr->indGr)) continue;
          for (int i = m_gr->indG[v]; i< m_gr->indG[v+1]; i++){
            int v2 = m_gr->nodeG[i];
            double w2 = m_gr->wG[i];


            if (Distances[v] + w2 < Distances[v2]) {
              Distances[v2] = Distances[v] + w2;
              if (add) distadd[v2] = distadd[v] + m_gr->add[i];
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
          if (add) distadd[i] = distadd[m_gr->nodeGr[j]] + m_gr->addr[j];
        }

      }
    }

    for (int i=0; i != m_arr.size(); i++){

      if (add) {
        m_result(k,i)= distadd[m_arr[i]];
      } else{
        m_result(k,i)= Distances[m_arr[i]];
      }


    }

    std::fill(Distances.begin(),Distances.end(),std::numeric_limits<double>::max());
    if (add){
      std::fill(distadd.begin(),distadd.end(),std::numeric_limits<double>::max());
    }

  }
}
