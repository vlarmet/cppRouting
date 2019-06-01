#include <iostream>
#include <queue>
#include <vector>
#include <fstream>
#include <limits>
#include <functional>
#include <Rcpp.h>
#include <cmath>



using namespace std;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::NumericVector Astar(std::vector<int> dep, std::vector<int> arr,std::vector<int> gfrom,std::vector<int> gto,std::vector<float> gw,int NbNodes,std::vector<float> lat,std::vector<float> lon,float k){
  
  
  Rcpp::NumericVector result(dep.size());
  
  
  
  
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
    int endNode=arr[j];
    float lata=lat[endNode];
    float lona=lon[endNode];
    
    std::vector<float> Distances(NbNodes, std::numeric_limits<float>::max()); 
    std::vector<float> Distances2(NbNodes, numeric_limits<float>::max());
    
    
    Distances[StartNode] = 0.0;                                                    
    Distances2[StartNode] = sqrt(pow(lat[StartNode]-lata,2)+pow(lon[StartNode]-lona,2))/k;
    
    std::vector<int> Parents(NbNodes, -1);                                             
   
    vector <int> closedList(NbNodes,0);
    vector <int> openList(NbNodes,0);
    priority_queue<std::pair<int, float>, vector<std::pair<int, float> >, comp > Q;
    Q.push(make_pair(StartNode,sqrt(pow(lat[StartNode]-lata,2)+pow(lon[StartNode]-lona,2))/k));                                             
    openList[StartNode]=1;
    
    while (!Q.empty()) {                                                          
      int v = Q.top().first;                                                      
      Q.pop();
      if (closedList[v]==1){
        continue;
      }
      openList[v]=0;
      closedList[v]=1;
      
      for (int i=0; i< G[v].size(); i++) {
        std::pair<int,float> j = G[v][i];                                                  
        int v2 = j.first;                                                     
        float w2 = j.second;
        if (closedList[v2]==1) {
          continue;
        }
        
        float temp;  
        temp = Distances[v] + w2;                              
        if (openList[v2]==0){
          
          Q.push(make_pair(v2,Distances2[v2]));
          openList[v2]=1;
        }
        
        
        else if (temp>=Distances[v2]){
          continue;
        }
        
        Parents[v2]=v;
        Distances[v2]=temp;
        Distances2[v2]=Distances[v2]+sqrt(pow(lat[v2]-lata,2)+pow(lon[v2]-lona,2))/k;
        Q.push(make_pair(v2,Distances2[v2]));
        openList[v2]=1;
      }
      
      
      if (v==endNode){
        break;
      }
      
    }
    
    if (Distances[endNode]==std::numeric_limits<float>::max()){
      result[j]= Rcpp::NumericVector::get_na();
    }
    else {
      result[j]= Distances[endNode];
    }
  
    
    
  }
  
  
  return result;
  
}


