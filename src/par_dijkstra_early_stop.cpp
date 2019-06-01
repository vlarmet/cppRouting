// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]

#include <iostream>
#include <queue>
#include <vector>
#include <fstream>
#include <limits>
#include <functional>
#include <Rcpp.h>
#include <string>
#include <RcppParallel.h>
#include <Rcpp.h>
using namespace RcppParallel;


struct Pardijkstra : public Worker
{
  //input
  const std::vector<std::vector<std::pair<int, float> > > m_graph;
  RVector<int> m_dep;
  RVector<int> m_arr;
  const int m_nbnodes;
  
  //output
  RcppParallel::RVector<double> m_result;
  
  //constructor
  Pardijkstra(const std::vector<std::vector<std::pair<int, float> > > graph,
              Rcpp::IntegerVector dep,
              Rcpp::IntegerVector arr,
              const int nbnodes,
              Rcpp::NumericVector result) : m_graph(graph),m_dep(dep), m_arr(arr),m_nbnodes(nbnodes),m_result(result)
  {
    
  }
  
  //overload () operator
  void operator()(std::size_t begin, std::size_t end){
    struct comp{
      
      bool operator()(const std::pair<int, float> &a, const std::pair<int, float> &b){
        return a.second > b.second;
      }
    };
    
    for (std::size_t k=begin; k!=end;k++){
      
      
      
      int StartNode=m_dep[k];
      
      std::vector<float> Distances(m_nbnodes, std::numeric_limits<float>::max());                   
      
      
      Distances[StartNode] = 0.0;                                                     
      
      std::vector<int> Parents(m_nbnodes, -1);                                            
      
      
      std::priority_queue<std::pair<int, float>, std::vector<std::pair<int, float> >, comp > Q;
      Q.push(std::make_pair(StartNode, 0.0));                                              
      
      while (!Q.empty()) {                                                        
        int v = Q.top().first;                                                     
        float w = Q.top().second;                                                     
        Q.pop();
        
        if (w <= Distances[v]) {                                                    
          
          for (int i=0; i< m_graph[v].size(); i++) {
            std::pair<int,float> j = m_graph[v][i];                                               
            int v2 = j.first;                                                      
            float w2 = j.second;
            
            if (Distances[v] + w2 < Distances[v2]) {                               
              Distances[v2] = Distances[v] + w2;                                   
              Parents[v2] = v;                                                   
              Q.push(std::make_pair(v2, Distances[v2]));
            }
          }
        }
        if (v==m_arr[k]){
          break;
        }
      }
      
      int EndNode=m_arr[k];
      
      if (Distances[EndNode]==std::numeric_limits<float>::max()){
        m_result[k] = Rcpp::NumericVector::get_na();
      }
      else {
        m_result[k] = Distances[EndNode];
      }
      
      
      
    }
    
  }
  
};


// [[Rcpp::export]]
Rcpp::NumericVector Dijkstra_early_stop_par(Rcpp::IntegerVector dep, Rcpp::IntegerVector arr,Rcpp::IntegerVector gfrom,Rcpp::IntegerVector gto,Rcpp::NumericVector gw,int NbNodes){
  
  
  Rcpp::NumericVector result(dep.size());
  int NbEdges=gfrom.size();
  
  std::vector<std::vector<std::pair<int, float> > > G(NbNodes);                                    
  for (int i = 0; i != NbEdges; ++i) {
    
    G[gfrom[i]].push_back(std::make_pair(gto[i], gw[i]));
    
    
  }
  
  
  Pardijkstra Dijfunc(G,dep,arr,NbNodes,result);
  parallelFor(0,dep.length(),Dijfunc);
  
  return Rcpp::wrap(result);
  
  
  
}


