// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]

#include <iostream>
#include <queue>
#include <vector>
#include <fstream>
#include <limits>
#include <functional>
#include <Rcpp.h>
#include <cmath>
#include <RcppParallel.h>
using namespace RcppParallel;


struct Pardijkstra : public Worker
{
  //input
  const std::vector<std::vector<std::pair<int, double> > > m_graph;
  RVector<int> m_dep;
  RVector<int> m_arr;
  const int m_nbnodes;
  RVector<double> m_lat;
  RVector<double> m_lon;
  const double m_k;
  
  //output
  RcppParallel::RVector<double> m_result;
  
  //constructor
  Pardijkstra(const std::vector<std::vector<std::pair<int, double> > > graph,
              Rcpp::IntegerVector dep,
              Rcpp::IntegerVector arr,
              const int nbnodes,
              Rcpp::NumericVector lat,
              Rcpp::NumericVector lon,
              const double k,
              Rcpp::NumericVector result) : m_graph(graph),m_dep(dep), m_arr(arr),m_nbnodes(nbnodes),m_lat(lat),m_lon(lon),m_k(k),m_result(result)
  {
    
  }
  
  //overload () operator
  void operator()(std::size_t begin, std::size_t end){
    struct comp{
      
      bool operator()(const std::pair<int, double> &a, const std::pair<int, double> &b){
        return a.second > b.second;
      }
    };
    
    for (std::size_t it=begin; it!=end;it++){
      
      
      
      int StartNode=m_dep[it];
      int EndNode=m_arr[it];
      double lata=m_lat[EndNode];
      double lona=m_lon[EndNode];
      
      std::vector<double> Distances(m_nbnodes, std::numeric_limits<double>::max());                   
      std::vector<double> Distances2(m_nbnodes, std::numeric_limits<double>::max());
      
      Distances[StartNode] = 0.0;
      Distances2[StartNode] = sqrt(pow(m_lat[StartNode]-lata,2)+pow(m_lon[StartNode]-lona,2))/m_k;
      
      std::vector<int> Parents(m_nbnodes, -1);                                            
      std::vector <int> closedList(m_nbnodes,0);
      std::vector <int> openList(m_nbnodes,0);
      
      std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double> >, comp > Q;
      Q.push(std::make_pair(StartNode,sqrt(pow(m_lat[StartNode]-lata,2)+pow(m_lon[StartNode]-lona,2))/m_k));                                           
      openList[StartNode]=1;
      
      while (!Q.empty()) {                                                    
        int v = Q.top().first;                                                
        Q.pop();
        if (closedList[v]==1){
          continue;
        }
        openList[v]=0;
        closedList[v]=1;
        
        for (int i=0; i< m_graph[v].size(); i++) {
          std::pair<int,double> j = m_graph[v][i];                                    
          int v2 = j.first;                                                      
          double w2 = j.second;
          if (closedList[v2]==1) {
            continue;
          }
          
          double temp;                             
          temp = Distances[v] + w2;
          if (openList[v2]==0){
            
            Q.push(std::make_pair(v2,Distances2[v2]));
            openList[v2]=1;
          }
          
          
          else if (temp>=Distances[v2]){
            continue;
          }
          
          Parents[v2]=v;
          Distances[v2]=temp;
          Distances2[v2]=Distances[v2]+sqrt(pow(m_lat[v2]-lata,2)+pow(m_lon[v2]-lona,2))/m_k;
          Q.push(std::make_pair(v2,Distances2[v2]));
          openList[v2]=1;
        }
        
        
        if (v==EndNode){
          break;
        }
        
      }
      
      if (Distances[EndNode]==std::numeric_limits<double>::max()){
        m_result[it] = Rcpp::NumericVector::get_na();
      }
      else {
        m_result[it] = Distances[EndNode];
      }
      
      
    }
    
  }
  
};


// [[Rcpp::export]]
Rcpp::NumericVector Astar_par(Rcpp::IntegerVector dep, Rcpp::IntegerVector arr,Rcpp::IntegerVector gfrom,Rcpp::IntegerVector gto,Rcpp::NumericVector gw,int NbNodes,Rcpp::NumericVector lat,Rcpp::NumericVector lon,double k){
  
  
  Rcpp::NumericVector result(dep.size());
  int NbEdges=gfrom.size();
  
  std::vector<std::vector<std::pair<int, double> > > G(NbNodes);                                    
  for (int i = 0; i != NbEdges; ++i) {
    
    G[gfrom[i]].push_back(std::make_pair(gto[i], gw[i]));
    
    
  }
  
  
  Pardijkstra Dijfunc(G,dep,arr,NbNodes,lat,lon,k,result);
  parallelFor(0,dep.length(),Dijfunc);
  
  return Rcpp::wrap(result);
  
  
  
}


