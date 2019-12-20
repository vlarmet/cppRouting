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

Rcpp::List Isochrone_multi(std::vector<int> gfrom,std::vector<int> gto,std::vector<double> gw,int NbNodes,std::vector<int> dep,std::vector<double> limit_vec,double max_limit,bool setdif,std::vector<std::string> dict,std::vector<int> keep){
  
  
  Rcpp::List finalresult(dep.size());
  
  
  
  
  struct comp{
    
    bool operator()(const std::pair<int, double> &a, const std::pair<int, double> &b){
      return a.second > b.second;
    }
  };
  
  
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
      if (Distances[v] > max_limit){
        break;
      }
    }
    
    std::vector<std::vector<std::string> > result(limit_vec.size());
    
    
    for (int i = 0; i < limit_vec.size(); i++){
      double lim=limit_vec[i];
      
      for (int k=0; k< Distances.size(); k++){
        if (keep[k]==1){
          if (Distances[k] < lim){
            if (setdif){
              Distances[k] = std::numeric_limits<double>::max();
            }
            
            result[i].push_back(dict[k]);
          }
        }
        

        
      
      }
    }
    
    finalresult[j] = result;
    //Reinitialize
    std::fill(Distances.begin(),Distances.end(),std::numeric_limits<double>::max());
    
    
    
  }
  
  
  return Rcpp::wrap(finalresult);
  
  
}


