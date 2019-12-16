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


struct ParBidir : public Worker
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
  RcppParallel::RVector<double> m_result;
  
  //constructor
  ParBidir(              const Rcpp::IntegerVector NodeG,
                            const Rcpp::NumericVector WG,
                            const Rcpp::IntegerVector IndG,
                            const Rcpp::IntegerVector NodeGr,
                            const Rcpp::NumericVector WGr,
                            const Rcpp::IntegerVector IndGr,
                            Rcpp::IntegerVector dep,
                            Rcpp::IntegerVector arr,
                            const int nbnodes,
                            Rcpp::NumericVector result) : m_NodeG(NodeG),m_WG(WG),m_IndG(IndG),
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
    std::vector<double> Distances2(m_nbnodes, std::numeric_limits<double>::max()); 
    std::vector <int> Visited(m_nbnodes,0);
    std::vector <int> Visiting(m_nbnodes,0);
    std::vector <int> Visited2(m_nbnodes,0);
    std::vector <int> Visiting2(m_nbnodes,0);
    for (std::size_t k=begin; k!=end;k++){
      
      
      
      int StartNode=m_dep[k];
      int EndNode=m_arr[k];
      

      
      
      Distances[StartNode] = 0.0;  
      Distances2[EndNode] = 0.0;

      
      std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double> >, comp > Q;
      std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double> >, comp > Qr;
      Q.push(std::make_pair(StartNode, 0.0)); 
      Qr.push(std::make_pair(EndNode, 0.0));
      Visiting[StartNode]=1;
      Visiting2[EndNode]=1;
      
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
            for (int i=m_IndG[v]; i< m_IndG[v+1]; i++){
              int v2 = m_NodeG[i];                                                     
              double w2 = m_WG[i];
              
              
              if (Distances[v] + w2 < Distances[v2]) {                               
                Distances[v2] = Distances[v] + w2;                                   
                                                               
                Q.push(std::make_pair(v2, Distances[v2]));
                Visiting[v2]=1;
                
              }
            }
          }
          if ((Visited2[v]==1 || Visiting2[v]==1)  && (Distances[v]+Distances2[v]) < mu){
            
            mu=Distances[v]+Distances2[v];
            
          }
        }
        
        if (!Qr.empty()){
          int vv=Qr.top().first;
          int ww=Qr.top().second;
          Qr.pop();
          
          
          Visited2[vv]=1;
          
          if (ww <= Distances2[vv]) {
            for (int i=m_IndGr[vv]; i< m_IndGr[vv+1]; i++){
              int vv2 = m_NodeGr[i];                                                     
              double ww2 = m_WGr[i];
              
              
              if (Distances2[vv] + ww2 < Distances2[vv2]) {                               
                Distances2[vv2] = Distances2[vv] + ww2;                                   
                                                               
                Qr.push(std::make_pair(vv2, Distances2[vv2]));
                Visiting2[vv2]=1;
              }
            }
          }
          if ((Visited[vv]==1 || Visiting[vv]==1) && (Distances[vv]+Distances2[vv]) < mu){
            
            mu=Distances[vv]+Distances2[vv];
            
          }
        }
        
      }
      
      
      
      if (mu >= std::numeric_limits<double>::max()){
        m_result[k] = Rcpp::NumericVector::get_na();
      }
      else {
        m_result[k]=mu;
      }
      
      std::fill(Distances.begin(), Distances.end(), std::numeric_limits<double>::max());
      std::fill(Distances2.begin(), Distances2.end(), std::numeric_limits<double>::max());
      std::fill(Visited.begin(), Visited.end(), 0);
      std::fill(Visited2.begin(), Visited2.end(), 0);
      std::fill(Visiting.begin(), Visiting.end(), 0);
      std::fill(Visiting2.begin(), Visiting2.end(), 0);
      
      
      
    }
    
  }
  
};


// [[Rcpp::export]]
Rcpp::NumericVector Bidir_par(Rcpp::IntegerVector dep, Rcpp::IntegerVector arr,Rcpp::IntegerVector gfrom,Rcpp::IntegerVector gto,Rcpp::NumericVector gw,int NbNodes){
  
  
  Rcpp::NumericVector result(dep.size());
  int NbEdges=gfrom.size();
  
  std::vector<std::vector<std::pair<int, double> > > G(NbNodes);   
  std::vector<std::vector<std::pair<int, double> > > Gr(NbNodes);
  int count=0;
  for (int i = 0; i != NbEdges; ++i) {
    
    G[gfrom[i]].push_back(std::make_pair(gto[i], gw[i]));
    Gr[gto[i]].push_back(std::make_pair(gfrom[i], gw[i]));
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
  
  
  ParBidir Dijfunc(NodeG,WG,IndG,NodeGr,WGr,IndGr,dep,arr,NbNodes,result);
  parallelFor(0,dep.length(),Dijfunc);
  
  return Rcpp::wrap(result);
  
  
  
}


