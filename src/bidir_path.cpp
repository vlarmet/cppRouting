#include <iostream>
#include <queue>
#include <vector>
#include <fstream>
#include <limits>
#include <functional>
#include <string>
#include <Rcpp.h>




using namespace std;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::List Bidir_path(std::vector<int> dep, std::vector<int> arr,std::vector<int> gfrom,std::vector<int> gto,std::vector<double> gw,int NbNodes,std::vector<std::string> dict,std::vector<int> keep){
  
  
  std::vector<std::vector<std::string> > result(dep.size());
  
  
  
  
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
  std::vector<int> Parents(NbNodes, -1);     
  std::vector<int> Parents2(NbNodes, -1);     
  vector <int> Visited(NbNodes,0);
  vector <int> Visiting(NbNodes,0);
  vector <int> Visited2(NbNodes,0);
  vector <int> Visiting2(NbNodes,0);
  for (unsigned int j=0; j!=dep.size();j++){
    if (j % 256){
      Rcpp::checkUserInterrupt ();
    }
    
    int StartNode=dep[j];
    int EndNode=arr[j];
    Distances[StartNode] = 0.0;  
    Distances2[EndNode] = 0.0;
    
    priority_queue<std::pair<int, double>, vector<std::pair<int, double> >, comp > Q;
    priority_queue<std::pair<int, double>, vector<std::pair<int, double> >, comp > Qr;
    Q.push(std::make_pair(StartNode, 0.0)); 
    Qr.push(std::make_pair(EndNode, 0.0)); 
    Visiting[StartNode]=1;
    Visiting2[EndNode]=1;
    
    int mid;
    double mu=std::numeric_limits<double>::max();
    
    while (!Q.empty() && !Qr.empty()) {  
      if (Q.top().second+Qr.top().second >= mu){
        break;
      }  
      
      
      if (!Q.empty()){
        int v=Q.top().first;
        int w=Q.top().second;
        Q.pop();
        
        
        Visited[v]=1;
        
        
        if (w <= Distances[v]) {
          for (int i=0; i< G[v].size(); i++){
            int v2 = G[v][i].first;                                                     
            double w2 = G[v][i].second;
            
            
            if (Distances[v] + w2 < Distances[v2]) {                               
              Distances[v2] = Distances[v] + w2;                                   
              Parents[v2] = v;                                                   
              Q.push(make_pair(v2, Distances[v2]));
              Visiting[v2]=1;
              
            }
          }
        }
        if ((Visited2[v]==1 || Visiting2[v]==1)  && (Distances[v]+Distances2[v]) < mu){
          mid=v;
          mu=Distances[v]+Distances2[v];
          
        }
      }
      
      if (!Qr.empty()){
        int vv=Qr.top().first;
        int ww=Qr.top().second;
        Qr.pop();
        
        
        Visited2[vv]=1;
        
        if (ww <= Distances2[vv]) {
          for (int i=0; i< Gr[vv].size(); i++){
            int vv2 = Gr[vv][i].first;                                                     
            double ww2 = Gr[vv][i].second;
            
            
            if (Distances2[vv] + ww2 < Distances2[vv2]) {                               
              Distances2[vv2] = Distances2[vv] + ww2;                                   
              Parents2[vv2] = vv;                                                   
              Qr.push(make_pair(vv2, Distances2[vv2]));
              Visiting2[vv2]=1;
            }
          }
        }
        if ((Visited[vv]==1 || Visiting[vv]==1) && (Distances[vv]+Distances2[vv]) < mu){
          mid=vv;
          mu=Distances[vv]+Distances2[vv];
          
        }
      }
      
    }
    
    
    
    std::vector <std::string> result2;
    for (auto p = Parents2[mid]; p != -1; p = Parents2[p]){
      if (keep[p]==1) result2.insert(result2.begin(),dict[p]);
    }
    
    if (Distances[mid]!=std::numeric_limits<double>::max() || Distances2[mid]!=std::numeric_limits<double>::max()){
      if (keep[mid]==1) result2.push_back(dict[mid]);
    }
    
    for (auto p = Parents[mid]; p != -1; p = Parents[p]){
      if (keep[p]==1) result2.push_back(dict[p]);
    }

    
 
    result[j] = result2;
    
    //Reinitialize
    std::fill(Distances.begin(),Distances.end(),std::numeric_limits<double>::max());
    std::fill(Distances2.begin(),Distances2.end(),std::numeric_limits<double>::max());
    std::fill(Visited.begin(),Visited.end(),0);
    std::fill(Visited2.begin(),Visited2.end(),0);
    std::fill(Visiting.begin(),Visiting.end(),0);
    std::fill(Visiting2.begin(),Visiting2.end(),0);
    std::fill(Parents.begin(),Parents.end(),-1);
    std::fill(Parents2.begin(),Parents2.end(),-1);
    
    
  }
  
  
  return Rcpp::wrap(result);
  
}


