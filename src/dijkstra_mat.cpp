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

Rcpp::NumericMatrix Dijkstra_mat(std::vector<int> gfrom,std::vector<int> gto,std::vector<double> gw,int NbNodes,std::vector<int> dep, std::vector<int> arr){
  
  
  //std::vector<std::vector<int> > result;
  Rcpp::NumericMatrix result(dep.size(),arr.size());
  
  
  
  struct comp{
    
    bool operator()(const std::pair<int, double> &a, const std::pair<int, double> &b){
      return a.second > b.second;
    }
  };
  
  //Graph
  
  int NbEdges=gfrom.size();
  
  std::vector<std::vector<std::pair<int, double> > > G(NbNodes);                                    
  
  for (int i = 0; i != NbEdges; ++i) {
    
    G[gfrom[i]].push_back(std::make_pair(gto[i], gw[i]));
    
    
  }
  
  //Boucle sur chaque trajet
  std::vector<double> Distances(NbNodes, std::numeric_limits<double>::max());
  for (int j=0; j!=dep.size();j++){
    if (j % 256){
      Rcpp::checkUserInterrupt ();
    }
    
    int StartNode=dep[j];
    
                       
   
    
    Distances[StartNode] = 0.0;                                                     
    
                                               
    
   
    priority_queue<std::pair<int, double>, vector<std::pair<int, double> >, comp > Q;
    Q.push(std::make_pair(StartNode, 0.0));                                             
    
    while (!Q.empty()) {                                                          
      int v = Q.top().first;                                                      
      double w = Q.top().second;                                                   
      Q.pop();
      
      if (w <= Distances[v]) {                                                    
        
        
        for (int i=0; i< G[v].size(); i++) {
          std::pair<int,double> j = G[v][i];                                                  
          int v2 = j.first;                                                      
          double w2 = j.second;
          
          if (Distances[v] + w2 < Distances[v2]) {                                
            Distances[v2] = Distances[v] + w2;                                    
                                                                
            Q.push(make_pair(v2, Distances[v2]));
          }
          
        }
        
        
      }
      
      
      
    }
    
    Rcpp::NumericVector result2(arr.size());
    for (int i=0;i!=arr.size();i++){
      if (Distances[arr[i]]==std::numeric_limits<double>::max()){
        result2[i] = Rcpp::NumericVector::get_na();
      }
      else {
        result2[i]= Distances[arr[i]];
      }
      
      
    }
    result.row(j) = result2;
    //Reinitialize
    std::fill(Distances.begin(),Distances.end(),std::numeric_limits<double>::max());
  }
  
  
  
  
  return result;
  
}


