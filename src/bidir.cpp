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

Rcpp::NumericVector Bidir(std::vector<int> dep, std::vector<int> arr,std::vector<int> gfrom,std::vector<int> gto,std::vector<double> gw,int NbNodes){
  
  
  Rcpp::NumericVector result(dep.size());
  
  
  
  
  struct comp{
    
    bool operator()(const std::pair<int, double> &a, const std::pair<int, double> &b){
      return a.second > b.second;
    }
  };
  
  //Graphs
  
  int NbEdges=gfrom.size();
  
  std::vector<std::vector<std::pair<int, double> > > G(NbNodes);   
  std::vector<std::vector<std::pair<int, double> > > Gr(NbNodes);
  
  for (unsigned int i = 0; i < NbEdges; ++i) {
    
    G[gfrom[i]].push_back(std::make_pair(gto[i], gw[i]));
    Gr[gto[i]].push_back(std::make_pair(gfrom[i], gw[i]));
    
    
  }
  
  
  //Boucle sur chaque trajet
  std::vector<double> Distances(NbNodes, std::numeric_limits<double>::max()); 
  std::vector<double> Distances2(NbNodes, std::numeric_limits<double>::max()); 
  std::vector <int> Visited(NbNodes,0);
  for (unsigned int k=0; k!=dep.size();k++){
    if (k % 256){
      Rcpp::checkUserInterrupt ();
    }
    
    int StartNode=dep[k];
    int EndNode=arr[k];
    
    Distances[StartNode] = 0.0;  
    Distances2[EndNode] = 0.0;
    std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double> >, comp > Q;
    std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double> >, comp > Qr;
    Q.push(std::make_pair(StartNode, 0.0)); 
    Qr.push(std::make_pair(EndNode, 0.0));
    Visited[StartNode]+=1;
    Visited[EndNode]+=1;
    
    double mu=std::numeric_limits<double>::max();
    
    
    while (!Q.empty() && !Qr.empty()) {  
      if (Q.top().second+Qr.top().second >= mu){
        break;
      }  
      
      
      if (!Q.empty()){
        int v=Q.top().first;
        int w=Q.top().second;
        Q.pop();
        
        
        //Visited[v]=1;
        
        
        if (w <= Distances[v]) {
          for (int i=0; i< G[v].size(); i++){
            int v2 = G[v][i].first;                                                     
            double w2 = G[v][i].second;
            
            
            if (Distances[v] + w2 < Distances[v2]) {                               
              Distances[v2] = Distances[v] + w2;                                   
              
              Q.push(std::make_pair(v2, Distances[v2]));
              Visited[v2]+=1;
              
            }
          }
        }
        if ((Visited[v]>1)  && (Distances[v]+Distances2[v]) < mu){
          
          mu=Distances[v]+Distances2[v];
          
        }
      }
      
      if (!Qr.empty()){
        int vv=Qr.top().first;
        int ww=Qr.top().second;
        Qr.pop();
        
        
        Visited[vv]+=1;
        
        if (ww <= Distances2[vv]) {
          for (int i=0; i< Gr[vv].size(); i++){
            int vv2 = Gr[vv][i].first;                                                     
            double ww2 = Gr[vv][i].second;
            
            
            if (Distances2[vv] + ww2 < Distances2[vv2]) {                               
              Distances2[vv2] = Distances2[vv] + ww2;                                   
              
              Qr.push(std::make_pair(vv2, Distances2[vv2]));
              Visited[vv2]+=1;
            }
          }
        }
        if ((Visited[vv]> 1) && (Distances[vv]+Distances2[vv]) < mu){
          
          mu=Distances[vv]+Distances2[vv];
          
        }
      }
      
    }
    
    
    
    if (mu >= std::numeric_limits<double>::max()){
      result[k] = Rcpp::NumericVector::get_na();
    }
    else {
      result[k]=mu;
    }
    
    //Reinitialize vectors
    std::fill(Distances.begin(),Distances.end(),std::numeric_limits<double>::max());
    std::fill(Distances2.begin(),Distances2.end(),std::numeric_limits<double>::max());
    std::fill(Visited.begin(),Visited.end(),0);

  }
  
  
  return result;
  
}


