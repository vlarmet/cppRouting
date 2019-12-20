// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]

#include <iostream>
#include <queue>
#include <vector>
#include <fstream>
#include <limits>
#include <functional>
#include <Rcpp.h>
#include <string>
#include <RcppParallel.h>
#include <Rcpp.h>
using namespace RcppParallel;


struct ParAstar : public Worker
{
  //input
  const RcppParallel::RVector<int> m_NodeG;
  const RcppParallel::RVector<double> m_WG;
  const RcppParallel::RVector<int> m_IndG;
  const RcppParallel::RVector<double> m_lat;
  const RcppParallel::RVector<double> m_lon;
  const double m_k;
  RVector<int> m_dep;
  RVector<int> m_arr;
  const int m_nbnodes;
  
  //output
  RcppParallel::RVector<double> m_result;
  
  //constructor
  ParAstar(              const Rcpp::IntegerVector NodeG,
                         const Rcpp::NumericVector WG,
                         const Rcpp::IntegerVector IndG,
                         const Rcpp::NumericVector lat,
                         const Rcpp::NumericVector lon,
                         const double k,
                         Rcpp::IntegerVector dep,
                         Rcpp::IntegerVector arr,
                         const int nbnodes,
                         Rcpp::NumericVector result) : m_NodeG(NodeG),m_WG(WG),m_IndG(IndG),
                         m_lat(lat),m_lon(lon),m_k(k),
                         m_dep(dep), m_arr(arr),m_nbnodes(nbnodes),m_result(result)
  {
    
  }
  
  //overload () operator
  void operator()(std::size_t begin, std::size_t end){
    struct comp{
      
      bool operator()(const std::pair<int, double> &a, const std::pair<int, double> &b){
        return a.second > b.second;
      }
    };
    
    
    std::vector<double> Distances(m_nbnodes, std::numeric_limits<double>::max()); 
    std::vector<double> Distances2(m_nbnodes, std::numeric_limits<double>::max()); 
    std::vector <int> closedList(m_nbnodes,0);
    std::vector <int> openList(m_nbnodes,0);
    for (std::size_t k=begin; k!=end;k++){
      
      int StartNode=m_dep[k];
      int endNode=m_arr[k];
      double lata=m_lat[endNode];
      double lona=m_lon[endNode];
      
      Distances[StartNode] = 0.0;                                                    
      Distances2[StartNode] = sqrt(pow(m_lat[StartNode]-lata,2)+pow(m_lon[StartNode]-lona,2))/m_k;
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
        
        for (int i=m_IndG[v]; i< m_IndG[v+1]; i++){
          int v2 = m_NodeG[i];                                                     
          double w2 = m_WG[i];
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
          
          
          Distances[v2]=temp;
          Distances2[v2]=Distances[v2]+sqrt(pow(m_lat[v2]-lata,2)+pow(m_lon[v2]-lona,2))/m_k;
          Q.push(std::make_pair(v2,Distances2[v2]));
          openList[v2]=1;
        }
        
        
        if (v==endNode){
          break;
        }
        
      }
      
      if (Distances[endNode]==std::numeric_limits<double>::max()){
        m_result[k]= Rcpp::NumericVector::get_na();
      }
      else {
        m_result[k]= Distances[endNode];
      }
      
      //Reinitialize vectors
      std::fill(Distances.begin(),Distances.end(),std::numeric_limits<double>::max());
      std::fill(Distances2.begin(),Distances2.end(),std::numeric_limits<double>::max());
      std::fill(closedList.begin(),closedList.end(),0);
      std::fill(openList.begin(),openList.end(),0);
      
      
    }
    
  }
  
};


// [[Rcpp::export]]
Rcpp::NumericVector Astar_par(Rcpp::IntegerVector dep, Rcpp::IntegerVector arr,Rcpp::IntegerVector gfrom,Rcpp::IntegerVector gto,Rcpp::NumericVector gw,int NbNodes,
                              Rcpp::NumericVector lat,Rcpp::NumericVector lon, double k){
  
  
  Rcpp::NumericVector result(dep.size());
  int NbEdges=gfrom.size();
  
  std::vector<std::vector<std::pair<int, double> > > G(NbNodes);   
  int count=0;
  for (int i = 0; i != NbEdges; ++i) {
    
    G[gfrom[i]].push_back(std::make_pair(gto[i], gw[i]));
    
    count+=1;
  }
  
  
  //Graph vectors
  Rcpp::IntegerVector NodeG(count);
  Rcpp::NumericVector WG(count);
  Rcpp::IntegerVector IndG(NbNodes+1);
  int count2=count;
  count=0;
  for (int i=0; i < G.size();i++){
    IndG[i]=count;
    
    for (int j=0; j < G[i].size();j++){
      NodeG[count]=G[i][j].first;
      WG[count]=G[i][j].second;
      count+=1;
    }
  }
  
  IndG[NbNodes]=count;
  std::vector<std::vector<std::pair<int, double> > > ().swap(G);

  
  
  ParAstar Dijfunc(NodeG,WG,IndG,lat,lon,k,dep,arr,NbNodes,result);
  parallelFor(0,dep.length(),Dijfunc);
  
  return Rcpp::wrap(result);
  
  
  
}


