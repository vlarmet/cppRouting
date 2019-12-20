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


struct ParPhast : public Worker
{
  //input
  const RcppParallel::RVector<int> m_NodeG;
  const RcppParallel::RVector<double> m_WG;
  const RcppParallel::RVector<int> m_IndG;
  const RcppParallel::RVector<int> m_NodeGr;
  const RcppParallel::RVector<double> m_WGr;
  const RcppParallel::RVector<int> m_IndGr;
  RVector<int> m_dep;
  RVector<int> m_arr;
  const int m_nbnodes;
  
  //output
  RcppParallel::RMatrix<double> m_result;
  
  //constructor
  ParPhast(              const Rcpp::IntegerVector NodeG,
                          const Rcpp::NumericVector WG,
                          const Rcpp::IntegerVector IndG,
                          const Rcpp::IntegerVector NodeGr,
                          const Rcpp::NumericVector WGr,
                          const Rcpp::IntegerVector IndGr,
                          Rcpp::IntegerVector dep,
                          Rcpp::IntegerVector arr,
                          const int nbnodes,
                          Rcpp::NumericMatrix result) : m_NodeG(NodeG),m_WG(WG),m_IndG(IndG),
                          m_NodeGr(NodeGr),m_WGr(WGr),m_IndGr(IndGr),
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
    
    for (std::size_t k=begin; k!=end;k++){
      
      int StartNode=m_dep[k];
      
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
          
          
          
          
          
          if (w <= Distances[v]) {
            for (int i=m_IndG[v]; i< m_IndG[v+1]; i++){
              int v2 = m_NodeG[i];                                                     
              double w2 = m_WG[i];
              
              
              if (Distances[v] + w2 < Distances[v2]) {                               
                Distances[v2] = Distances[v] + w2;                                   
                Q.push(std::make_pair(v2, Distances[v2]));
                
                
              }
            }
          }
          
        }
        
        
      }
      
      
      
      //Backward
      
      for (int i=0; i < (m_IndGr.size()-1); i++){
        
        for (int j=m_IndGr[i]; j < m_IndGr[i+1]; j++){
          
          
          if (Distances[m_NodeGr[j]]+m_WGr[j] < Distances[i]) Distances[i] = Distances[m_NodeGr[j]]+m_WGr[j];
          
          
          
        }
        
        
      }
      
      
      for (int i=0;i!=m_arr.size();i++){
        if (Distances[m_arr[i]]==std::numeric_limits<double>::max()){
          m_result(k,i) = Rcpp::NumericVector::get_na();
        }
        else {
          m_result(k,i)= Distances[m_arr[i]];
        }
        
        
      }
      
      
      std::fill(Distances.begin(),Distances.end(),std::numeric_limits<double>::max());
      
    }
      
      
    }
    
  
  
};


// [[Rcpp::export]]
Rcpp::NumericMatrix Phast_par(Rcpp::IntegerVector dep, Rcpp::IntegerVector arr,Rcpp::IntegerVector gfrom,Rcpp::IntegerVector gto,Rcpp::NumericVector gw,int NbNodes){
  
  Rcpp::NumericMatrix result(dep.size(),arr.size());
  int NbEdges=gfrom.size();
  
  std::vector<std::vector<std::pair<int, double> > > G(NbNodes);   
  std::vector<std::vector<std::pair<int, double> > > Gr(NbNodes);
  int count=0;
  for (unsigned int i = 0; i < NbEdges; ++i) {
    
    if (gfrom[i] > gto[i])  G[gfrom[i]].push_back(std::make_pair(gto[i], gw[i]));
    if (gfrom[i] < gto[i]) Gr[gto[i]].push_back(std::make_pair(gfrom[i], gw[i]));
    count+=1;
    
  }
  for (int i=0; i < Gr.size(); i++){
    std::sort(Gr[i].begin(),Gr[i].end());
    //std::reverse(Gr[i].begin(),Gr[i].end());
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
  //Buckets
  
  Rcpp::IntegerVector NodeGr(count2);
  Rcpp::NumericVector WGr(count2);
  Rcpp::IntegerVector IndGr(NbNodes+1);
  count=0;
  for (int i=0; i < Gr.size();i++){
    IndGr[i]=count;
    
    for (int j=0; j < Gr[i].size();j++){
      NodeGr[count]=Gr[i][j].first;
      WGr[count]=Gr[i][j].second;
      count+=1;
    }
  }
  IndGr[NbNodes]=count;
  std::vector<std::vector<std::pair<int, double> > > ().swap(Gr);
  
  
  ParPhast Dijfunc(NodeG,WG,IndG,NodeGr,WGr,IndGr,dep,arr,NbNodes,result);
  
  parallelFor(0,dep.length(),Dijfunc);
  
  return Rcpp::wrap(result);
  
  
  
}


