#include <iostream>
#include <queue>
#include <vector>
#include <fstream>
#include <limits>
#include <functional>
#include <chrono>
#include <thread>
#include <Rcpp.h>




using namespace std;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::NumericVector Bidir_mod3(std::vector<int> &dep, std::vector<int> &arr,std::vector<int> &gfrom,std::vector<int> &gto,std::vector<double> &gw,int NbNodes,std::vector<int> &Rank){
  
  
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
    
    
   if (Rank[gfrom[i]] < Rank[gto[i]])  G[gfrom[i]].push_back(std::make_pair(gto[i], gw[i]));
   if (Rank[gfrom[i]] > Rank[gto[i]]) Gr[gto[i]].push_back(std::make_pair(gfrom[i], gw[i]));
    
    
  }
  
  
  //Boucle sur chaque trajet
  
  std::vector<double> Distances(NbNodes, std::numeric_limits<double>::max()); 
  std::vector<double> Distances2(NbNodes, std::numeric_limits<double>::max()); 
  std::vector <int> Visited1(NbNodes,0);
  std::vector <int> Visited2(NbNodes,0);
  std::vector <int> Visited;
  
  
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
    
    double mu=std::numeric_limits<double>::max();
    
    
    
    //Rcpp::Rcerr <<k<<std::endl;
    
    
    
    while (true) {  
      
      
      
      if (Q.top().second > mu && Qr.top().second > mu){
        break;
      }  
      if (Q.empty() && Qr.empty()){
        break;
      }  
      
      
      
      if (!Q.empty()){
        
        int v=Q.top().first;
        double w=Q.top().second;
        Q.pop();
        
        
        
        
        Visited.push_back(v);
        
        Visited1[v]=1;
        
        if ((Visited2[v]==1)  && (Distances[v]+Distances2[v]) < mu){
          
          mu=Distances[v]+Distances2[v];
          
        }
        
        
        
        
        
        if (w <= Distances[v]) {
          for (int i=0; i< G[v].size(); i++){
            int v2 = G[v][i].first;                                                     
            double w2 = G[v][i].second;
            
            //if (Rank[v] > Rank[v2]) continue;
            
            if (Distances[v] + w2 < Distances[v2]) {                               
              Distances[v2] = Distances[v] + w2;                                   
              Q.push(std::make_pair(v2, Distances[v2]));
              Visited1[v2]=1;
              
              Visited.push_back(v2);
              
            }
          }
        }
        
      }
      
      if (!Qr.empty()){
        int vv=Qr.top().first;
        double ww=Qr.top().second;
        Qr.pop();
        
        Visited.push_back(vv);
        
        Visited2[vv]=1;
        
        
        if ((Visited1[vv]== 1) && (Distances[vv]+Distances2[vv]) < mu){
          
          mu=Distances[vv]+Distances2[vv];
          
        }
        
        if (ww <= Distances2[vv]) {
          for (int i=0; i< Gr[vv].size(); i++){
            int vv2 = Gr[vv][i].first;                                                     
            double ww2 = Gr[vv][i].second;
            
            //if (Rank[vv] > Rank[vv2]) continue;
            
            if (Distances2[vv] + ww2 < Distances2[vv2]) {                               
              Distances2[vv2] = Distances2[vv] + ww2;                                   
              
              Qr.push(std::make_pair(vv2, Distances2[vv2]));
              Visited2[vv2]=1;
              
              Visited.push_back(vv2);
              
            }
          }
        }
        
      }
      
    }

    for (int i=0; i<Visited.size();i++){
      Distances[Visited[i]]=std::numeric_limits<double>::max();
      Distances2[Visited[i]]=std::numeric_limits<double>::max();
      Visited1[Visited[i]]=0;
      Visited2[Visited[i]]=0;
    }
    Visited.clear();
    
    
    if (mu >= std::numeric_limits<double>::max()){
      result[k] = Rcpp::NumericVector::get_na();
    }
    else {
      result[k]=mu;
      
    }
    
    
  }
  
  
  return result;
  
}


