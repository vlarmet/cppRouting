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

struct ParBidirmod : public Worker
{
  //input
  //const std::vector<std::vector<std::pair<int,double> > > &m_Bucket;
  //const std::vector<std::vector<std::pair<int,double> > > &m_Graphr;
  const RcppParallel::RVector<int> m_NodeG;
  const RcppParallel::RVector<double> m_WG;
  const RcppParallel::RVector<int> m_IndG;
  const RcppParallel::RVector<int> m_NodeGr;
  const RcppParallel::RVector<double> m_WGr;
  const RcppParallel::RVector<int> m_IndGr;
  
  const RcppParallel::RVector<int> m_Dep;
  const RcppParallel::RVector<int> m_Arr;
  const int m_Nb;
  
  
  //output
  RcppParallel::RVector<double> m_result;
  
  //constructor
  ParBidirmod(//const std::vector<std::vector<std::pair<int,double> > > &Bucket,
    //const std::vector<std::vector<std::pair<int,double> > > &Graphr,
    const Rcpp::IntegerVector NodeG,
    const Rcpp::NumericVector WG,
    const Rcpp::IntegerVector IndG,
    const Rcpp::IntegerVector NodeGr,
    const Rcpp::NumericVector WGr,
    const Rcpp::IntegerVector IndGr,
    const Rcpp::IntegerVector Dep,
    const Rcpp::IntegerVector Arr,
    const int Nb,
    Rcpp::NumericVector Result) : m_NodeG(NodeG),m_WG(WG),m_IndG(IndG),
    m_NodeGr(NodeGr),m_WGr(WGr),m_IndGr(IndGr),
    m_Dep(Dep),m_Arr(Arr),m_Nb(Nb),m_result(Result)
  {
    
  }
  
  //overload () operator
  void operator()(std::size_t begin, std::size_t end){
    struct comp{
      
      bool operator()(const std::pair<int, double> &a, const std::pair<int, double> &b){
        return a.second > b.second;
      }
    };
    
    std::vector<double> Distances(m_Nb, std::numeric_limits<double>::max()); 
    std::vector<double> Distances2(m_Nb, std::numeric_limits<double>::max()); 
    std::vector<int> Visited1(m_Nb,0);
    std::vector<int> Visited2(m_Nb,0);
    std::vector<int> Visited;
    
    double mu=std::numeric_limits<double>::max();
    
    
    
    for (std::size_t k=begin; k!=end;k++){
      
      unsigned int StartNode = static_cast <unsigned int> (m_Dep [k]);
      unsigned int EndNode = static_cast <unsigned int> (m_Arr [k]);
      //int StartNode=m_Dep[k];
      //int EndNode=m_Arr[k];
      
      Distances[StartNode] = 0.0;  
      Distances2[EndNode] = 0.0;
      
      std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double> >, comp > Q;
      std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double> >, comp > Qr;
      Q.push(std::make_pair(StartNode, 0.0)); 
      Qr.push(std::make_pair(EndNode, 0.0));
      
      
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
            for (int i=m_IndG[v]; i< m_IndG[v+1]; i++){
              int v2 = m_NodeG[i];                                                     
              double w2 = m_WG[i];
              
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
            for (int i=m_IndGr[vv]; i< m_IndGr[vv+1]; i++){
              int vv2 = m_NodeGr[i];                                                     
              double ww2 = m_WGr[i];
              
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
      
      
      if (mu >= std::numeric_limits<double>::max()){
        m_result[k] = Rcpp::NumericVector::get_na();
      }
      else {
        m_result[k]=mu;
        
      }
      for (int i=0; i<Visited.size();i++) Distances[Visited[i]]=std::numeric_limits<double>::max();
      
      for (int i=0; i<Visited.size();i++) Distances2[Visited[i]]=std::numeric_limits<double>::max();
      for (int i=0; i<Visited.size();i++) Visited1[Visited[i]]=0;
      for (int i=0; i<Visited.size();i++) Visited2[Visited[i]]=0;
      
      Visited.clear();
      
      mu=std::numeric_limits<double>::max();
      
      
      
      
    }
    
    
  }
  
};





// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::NumericVector par_Bidir_mod2(std::vector<int> &dep, std::vector<int> &arr,std::vector<int> &gfrom,std::vector<int> &gto,std::vector<double> &gw,int NbNodes,std::vector<int> &Rank){
  
  
  struct comp{
    
    bool operator()(const std::pair<int, double> &a, const std::pair<int, double> &b){
      return a.second > b.second;
    }
  };
  
  //Graphs
  
  int NbEdges=gfrom.size();
  
  std::vector<std::vector<std::pair<int, double> > > G(NbNodes);   
  std::vector<std::vector<std::pair<int, double> > > Gr(NbNodes);
  
  int count=0;
  int count2=0;
  for (unsigned int i = 0; i < NbEdges; ++i) {
    
    if (Rank[gfrom[i]] < Rank[gto[i]]){
      G[gfrom[i]].push_back(std::make_pair(gto[i], gw[i]));
      count+=1;
    }
    if (Rank[gfrom[i]] > Rank[gto[i]]) {
      Gr[gto[i]].push_back(std::make_pair(gfrom[i], gw[i]));
      count2+=1;
    }
    
    
    
  }
  
  
  //Graph vectors
  Rcpp::IntegerVector NodeG(count);
  Rcpp::NumericVector WG(count);
  Rcpp::IntegerVector IndG(NbNodes+1);
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
  
  Rcpp::NumericVector result(arr.size());
  Rcpp::IntegerVector dep2=Rcpp::wrap(dep);
  Rcpp::IntegerVector arr2=Rcpp::wrap(arr);
  
  ParBidirmod Dijfunc(NodeG,WG,IndG,NodeGr,WGr,IndGr,dep2,arr2,NbNodes,result);
  parallelFor(0,arr2.length(),Dijfunc);
  
  
  return result;
  
  
}


