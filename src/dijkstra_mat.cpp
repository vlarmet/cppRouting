#include <iostream>
#include <queue>
#include <vector>
#include <fstream>
#include <limits>
#include <functional>
#include <Rcpp.h>
#include <string>



using namespace std;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::NumericMatrix Dijkstra_mat(std::vector<int> gfrom,std::vector<int> gto,std::vector<float> gw,int NbNodes,std::vector<int> dep, std::vector<int> arr){
  
  
  //std::vector<std::vector<int> > result;
  Rcpp::NumericMatrix result(dep.size(),arr.size());
  
  
  
  struct comp{
    
    bool operator()(const std::pair<int, float> &a, const std::pair<int, float> &b){
      return a.second > b.second;
    }
  };
  
  //Graph
  
  int NbEdges=gfrom.size();
  
  std::vector<std::vector<std::pair<int, float> > > G(NbNodes);                                    
  
  for (int i = 0; i != NbEdges; ++i) {
    
    G[gfrom[i]].push_back(std::make_pair(gto[i], gw[i]));
    
    
  }
  
  //Boucle sur chaque trajet
  
  for (int j=0; j!=dep.size();j++){
    if (j % 256){
      Rcpp::checkUserInterrupt ();
    }
    
    int StartNode=dep[j];
    
    std::vector<float> Distances(NbNodes, std::numeric_limits<float>::max());                   
   
    
    Distances[StartNode] = 0;                                                     
    
    std::vector<int> Parents(NbNodes, -1);                                             
    
   
    priority_queue<std::pair<int, float>, vector<std::pair<int, float> >, comp > Q;
    Q.push(std::make_pair(StartNode, 0));                                             
    
    while (!Q.empty()) {                                                          
      int v = Q.top().first;                                                      
      float w = Q.top().second;                                                   
      Q.pop();
      
      if (w <= Distances[v]) {                                                    
        
        
        for (int i=0; i< G[v].size(); i++) {
          std::pair<int,float> j = G[v][i];                                                  
          int v2 = j.first;                                                      
          float w2 = j.second;
          
          if (Distances[v] + w2 < Distances[v2]) {                                
            Distances[v2] = Distances[v] + w2;                                    
            Parents[v2] = v;                                                      
            Q.push(make_pair(v2, Distances[v2]));
          }
          
        }
        
        
      }
      
      
      
    }
    
    Rcpp::NumericVector result2(arr.size());
    for (int i=0;i!=arr.size();i++){
      if (Distances[arr[i]]==std::numeric_limits<float>::max()){
        result2[i] = Rcpp::NumericVector::get_na();
      }
      else {
        result2[i]= Distances[arr[i]];
      }
      
      
    }
    result.row(j) = result2;
    
  }
  
  
  
  
  return result;
  
}


