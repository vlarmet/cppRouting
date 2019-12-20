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

Rcpp::NumericMatrix Bidir_mat3(std::vector<int> &dep, std::vector<int> &arr,std::vector<int> &gfrom,std::vector<int> &gto,std::vector<double> &gw,int NbNodes,std::vector<int> &Rank){
  
  Rcpp::NumericMatrix result(arr.size(),dep.size());

  
  std::vector<std::vector<std::pair<int,double> > > BucketF(NbNodes);
  
  
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
  
  
  //Forward
  
  std::vector<double> Distances(NbNodes, std::numeric_limits<double>::max()); 
  std::vector <int> Visited(NbNodes,0);

  
  
  for (unsigned int k=0; k!=dep.size();k++){
    if (k % 256){
      Rcpp::checkUserInterrupt ();
    }
    
    auto t1 = std::chrono::high_resolution_clock::now();
    int StartNode=dep[k];

    Distances[StartNode] = 0.0;  
    
    
    std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double> >, comp > Q;
    Q.push(std::make_pair(StartNode, 0.0)); 

    
    while (true) {  
      
      if (Q.empty()){
        break;
      }  
      
      
      
      if (!Q.empty()){
        
        int v=Q.top().first;
        double w=Q.top().second;
        Q.pop();

        
        Visited[v]=1;
        
        
        if (w <= Distances[v]) {
          for (int i=0; i< G[v].size(); i++){
            int v2 = G[v][i].first;                                                     
            double w2 = G[v][i].second;

            
            if (Distances[v] + w2 < Distances[v2]) {                               
              Distances[v2] = Distances[v] + w2;                                   
              Q.push(std::make_pair(v2, Distances[v2]));
              Visited[v2]=1;

            }
          }
        }
        
      }
      
      
    }
    
    for (int i=0; i < Distances.size();i++){
      if (Distances[i]<std::numeric_limits<double>::max()){
        BucketF[i].push_back(std::make_pair(k,Distances[i]));
      }
      
    }
    

    
    std::fill(Distances.begin(), Distances.end(), std::numeric_limits<double>::max());

  }
  for (int i=0; i < BucketF.size();i++){
    std::sort(BucketF[i].begin(),BucketF[i].end());
  }

  //Backward
  
  
  for (unsigned int k=0; k!=arr.size();k++){
    if (k % 256){
      Rcpp::checkUserInterrupt ();
    }
    

    
    
    Rcpp::NumericVector Dist(dep.size(), std::numeric_limits<double>::max());
    for (int i=0; i < dep.size();i++) Dist[i]=std::numeric_limits<double>::max();
    
    
    
    int StartNode=arr[k];

    Distances[StartNode] = 0.0;  
    
    
    std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double> >, comp > Q;
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
        
        if (Visited[v]>0){
          
          for (int i=0; i < BucketF[v].size();i++){
            int Source=BucketF[v][i].first;
            
          
            double D=Distances[v] + BucketF[v][i].second;
            
            
            
            if (Dist[Source] > D) {
              
              Dist[Source] = D;
            }
   
        }
        }

        if (w <= Distances[v]) {
          for (int i=0; i< Gr[v].size(); i++){
            int v2 = Gr[v][i].first;                                                     
            double w2 = Gr[v][i].second;
            
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
      if (Dist[i]==std::numeric_limits<double>::max()) Dist[i]=Rcpp::NumericVector::get_na();
    }
    result.row(k) = Dist;
    
    
    
    
  }

  
  return Rcpp::transpose(result);

  
}


