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

Rcpp::NumericVector Astar(std::vector<int> dep, std::vector<int> arr,std::vector<int> gfrom,std::vector<int> gto,std::vector<double> gw,int NbNodes,std::vector<double> lat,std::vector<double> lon,double k){
  
  
  Rcpp::NumericVector result(dep.size());
  
  
  
  
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
  std::vector<double> Distances2(NbNodes, numeric_limits<double>::max());
  vector <int> closedList(NbNodes,0);
  vector <int> openList(NbNodes,0);
  
  for (int j=0; j!=dep.size();j++){
    if (j % 256){
      Rcpp::checkUserInterrupt ();
    }
    
    int StartNode=dep[j];
    int endNode=arr[j];
    double lata=lat[endNode];
    double lona=lon[endNode];
    
    Distances[StartNode] = 0.0;                                                    
    Distances2[StartNode] = sqrt(pow(lat[StartNode]-lata,2)+pow(lon[StartNode]-lona,2))/k;
    priority_queue<std::pair<int, double>, vector<std::pair<int, double> >, comp > Q;
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
        std::pair<int,double> j = G[v][i];                                                  
        int v2 = j.first;                                                     
        double w2 = j.second;
        if (closedList[v2]==1) {
          continue;
        }
        
        double temp;  
        temp = Distances[v] + w2;                              
        if (openList[v2]==0){
          
          Q.push(make_pair(v2,Distances2[v2]));
          openList[v2]=1;
        }
        
        
        else if (temp>=Distances[v2]){
          continue;
        }
        
        
        Distances[v2]=temp;
        Distances2[v2]=Distances[v2]+sqrt(pow(lat[v2]-lata,2)+pow(lon[v2]-lona,2))/k;
        Q.push(make_pair(v2,Distances2[v2]));
        openList[v2]=1;
      }
      
      
      if (v==endNode){
        break;
      }
      
    }
    
    if (Distances[endNode]==std::numeric_limits<double>::max()){
      result[j]= Rcpp::NumericVector::get_na();
    }
    else {
      result[j]= Distances[endNode];
    }
  
  //Reinitialize vectors
  std::fill(Distances.begin(),Distances.end(),std::numeric_limits<double>::max());
  std::fill(Distances2.begin(),Distances2.end(),std::numeric_limits<double>::max());
  std::fill(closedList.begin(),closedList.end(),0);
  std::fill(openList.begin(),openList.end(),0);
    
    
  }
  
  
  return result;
  
}


