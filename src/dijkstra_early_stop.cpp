#include <iostream>
#include <queue>
#include <vector>
#include <fstream>
#include <limits>
#include <functional>
#include <Rcpp.h>




using namespace std;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::NumericVector Dijkstra_early_stop(std::vector<int> dep, std::vector<int> arr,std::vector<int> gfrom,std::vector<int> gto,std::vector<double> gw,int NbNodes){
  
  
  Rcpp::NumericVector result(dep.size());
  
  
  
  
  struct comp{
    
    bool operator()(const std::pair<int, double> &a, const std::pair<int, double> &b){
      return a.second > b.second;
    }
  };
  
  //Graph
  
  int NbEdges=gfrom.size();
  
  std::vector<std::vector<std::pair<int, double> > > G(NbNodes);                                    
  
  for (unsigned int i = 0; i < NbEdges; ++i) {
    
    G[gfrom[i]].push_back(std::make_pair(gto[i], gw[i]));
    
    
  }
  //Rcpp::Rcerr << "Graph construit!\n";
  
  //Boucle sur chaque trajet
  std::vector<double> Distances(NbNodes, std::numeric_limits<double>::max());      
  
  for (unsigned int j=0; j!=dep.size();j++){
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
        
        for (unsigned int i=0; i< G[v].size(); i++) {
          std::pair<int,double> j = G[v][i];                                               
          int v2 = j.first;                                                      
          double w2 = j.second;
          
          if (Distances[v] + w2 < Distances[v2]) {                               
            Distances[v2] = Distances[v] + w2;                                   
                                                              
            Q.push(make_pair(v2, Distances[v2]));
          }
        }
      }
      if (v==arr[j]){
        break;
      }
    }
    
    int EndNode=arr[j];
    if (Distances[EndNode]==std::numeric_limits<double>::max()){
      result[j] = Rcpp::NumericVector::get_na();
    }
    else {
      result[j]=Distances[EndNode];
    }
    
    //Reinitialize
    std::fill(Distances.begin(),Distances.end(),std::numeric_limits<double>::max());
    
  }
  
  
  return result;
  
}


