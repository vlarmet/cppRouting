#include <iostream>
#include <queue>
#include <vector>
#include <fstream>
#include <limits>
#include <string>
#include <functional>
#include <Rcpp.h>




using namespace std;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::List NBA_path(std::vector<int> dep, std::vector<int> arr,std::vector<int> gfrom,std::vector<int> gto,std::vector<double> gw,int NbNodes,std::vector<double> lat,std::vector<double> lon,double k,std::vector<std::string> dict,std::vector<int> keep){
  
  
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
  vector <int> Visited(NbNodes,0);
  vector <int> Visited1check(NbNodes,0);
  vector <int> Visited2check(NbNodes,0);
  std::vector<int> Parents(NbNodes, -1);     
  std::vector<int> Parents2(NbNodes, -1); 
  for (unsigned int j=0; j!=dep.size();j++){
    if (j % 256){
      Rcpp::checkUserInterrupt ();
    }
    
    int StartNode=dep[j];
    int EndNode=arr[j];
    double lata=lat[EndNode];
    double lona=lon[EndNode];
    double lata2=lat[StartNode];
    double lona2=lon[StartNode];
    
    Distances[StartNode] = 0.0;  
    Visited1check[StartNode]=1;
    Distances2[EndNode] = 0.0;
    Visited2check[EndNode]=1;
    

    priority_queue<std::pair<int, double>, vector<std::pair<int, double> >, comp > Q;
    priority_queue<std::pair<int, double>, vector<std::pair<int, double> >, comp > Qr;
    Q.push(std::make_pair(StartNode, sqrt(pow(lat[StartNode]-lata,2)+pow(lon[StartNode]-lona,2))/k)); 
    Qr.push(std::make_pair(EndNode, sqrt(pow(lat[EndNode]-lata2,2)+pow(lon[EndNode]-lona2,2))/k)); 
    
    
    
    //double Pr=0.5*sqrt(pow(lat[StartNode]-lata,2)+pow(lon[StartNode]-lona,2))/k;
    double total1=sqrt(pow(lata2-lata,2)+pow(lona2-lona,2))/k;
    double total2=sqrt(pow(lata2-lata,2)+pow(lona2-lona,2))/k;
    int mid;
    double mu=std::numeric_limits<double>::max();
    
    while (!Q.empty() && !Qr.empty()) {  
      //Forward
      if (Q.size() < Qr.size()){
        int v=Q.top().first;
        Q.pop();
        if (Visited[v]==0){
          Visited[v]=1;
          
          if ((Distances[v] + sqrt(pow(lat[v]-lata,2)+pow(lon[v]-lona,2))/k) >= mu || (Distances[v] + total2 - sqrt(pow(lat[v]-lata2,2)+pow(lon[v]-lona2,2))/k) >= mu){}
          
          else {
            for (int i=0; i < G[v].size(); i++){
              int v2=G[v][i].first;
              double w2=G[v][i].second;
              if (Visited[v2]==1){
                continue;
              }
              double tentative=Distances[v]+w2;
              
              if (Visited1check[v2]==0  || Distances[v2] > tentative){
                Distances[v2]=tentative;
                Visited1check[v2]=1;
                Parents[v2]=v;
                Q.push(std::make_pair(v2, tentative + sqrt(pow(lat[v2]-lata,2)+pow(lon[v2]-lona,2))/k));
                
                if (Visited2check[v2]==1){
                  double temp=tentative + Distances2[v2];
                  if (mu > temp){
                    mu=temp;
                    mid=v2;
                    
                  }
                }
                
              }
              
            }
          }
          
          if (!Q.empty()){
            total1=Q.top().second;
          }
          
          
        } 
      }
      
      //Backward
      else {
        
        int vv=Qr.top().first;
        
        Qr.pop();
        
        if (Visited[vv]==0){
          
          Visited[vv]=1;
          
          if ((Distances2[vv] + sqrt(pow(lat[vv]-lata2,2)+pow(lon[vv]-lona2,2))/k) >= mu || (Distances2[vv] + total1 - sqrt(pow(lat[vv]-lata,2)+pow(lon[vv]-lona,2))/k) >= mu){}
          
          else {
            for (int i=0; i < Gr[vv].size(); i++){
              int vv2=Gr[vv][i].first;
              double ww2=Gr[vv][i].second;
              if (Visited[vv2]==1){
                continue;
              }
              double tentative=Distances2[vv]+ww2;
              
              if (Visited2check[vv2]==0  || Distances2[vv2] > tentative){
                Distances2[vv2]=tentative;
                Visited2check[vv2]=1;
                Parents2[vv2]=vv;
                Qr.push(std::make_pair(vv2, tentative + sqrt(pow(lat[vv2]-lata2,2)+pow(lon[vv2]-lona2,2))/k));
                
                if (Visited1check[vv2]==1){
                  double temp=tentative + Distances[vv2];
                  if (mu > temp){
                    mu=temp;
                    mid=vv2;
                    
                  }
                }
                
              }
              
            }
          }
          
          if (!Qr.empty()){
            total2=Qr.top().second;
          }
          
        } 
      }
      
      
    }
    
    
    std::vector <std::string> result2;
    for (auto p = Parents2[mid]; p != -1; p = Parents2[p]){
      
    if (keep[p]==1)  result2.insert(result2.begin(),dict[p]);
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
    std::fill(Visited1check.begin(),Visited1check.end(),0);
    std::fill(Visited2check.begin(),Visited2check.end(),0);
    std::fill(Parents.begin(),Parents.end(),-1);
    std::fill(Parents2.begin(),Parents2.end(),-1);
    
    
  }
  
  
  return Rcpp::wrap(result);
  
}


