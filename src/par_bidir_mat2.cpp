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

struct ParBidirMat : public Worker
{
  //input
  //const std::vector<std::vector<std::pair<int,double> > > &m_Bucket;
  //const std::vector<std::vector<std::pair<int,double> > > &m_Graphr;
  const RcppParallel::RVector<int> m_NodeGr;
  const RcppParallel::RVector<double> m_WGr;
  const RcppParallel::RVector<int> m_IndGr;
  const RcppParallel::RVector<int> m_NodeBu;
  const RcppParallel::RVector<double> m_WBu;
  const RcppParallel::RVector<int> m_IndBu;
  
  const RcppParallel::RVector<int> m_Arr;
  const int m_Size;
  const int m_Nb;
  
  
  //output
  RcppParallel::RMatrix<double> m_result;
  
  //constructor
  ParBidirMat(//const std::vector<std::vector<std::pair<int,double> > > &Bucket,
              //const std::vector<std::vector<std::pair<int,double> > > &Graphr,
              const Rcpp::IntegerVector NodeGr,
              const Rcpp::NumericVector WGr,
              const Rcpp::IntegerVector IndGr,
              const Rcpp::IntegerVector NodeBu,
              const Rcpp::NumericVector WBu,
              const Rcpp::IntegerVector IndBu,
              const Rcpp::IntegerVector Arr,
              const int Size,
              const int Nb,
              Rcpp::NumericMatrix Result) : m_NodeGr(NodeGr),m_WGr(WGr),m_IndGr(IndGr),
              m_NodeBu(NodeBu),m_WBu(WBu),m_IndBu(IndBu),
              m_Arr(Arr),m_Size(Size),m_Nb(Nb),m_result(Result)
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
    for (std::size_t k=begin; k!=end;k++){
      
      //Rcpp::Rcerr << k <<std::endl;
      
      std::vector<double> Dist(m_Size,std::numeric_limits<double>::max());
      
      int StartNode=m_Arr[k];
      
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
          
          if ((m_IndBu[v+1] - m_IndBu[v]) > 0){
            
            for (int i=m_IndBu[v]; i < m_IndBu[v+1];i++){
              
              
              int Source=m_NodeBu[i];
              
              
              double D=Distances[v] + m_WBu[i];
              
              if (Dist[Source] > D) {
                
                Dist[Source] = D;
              }
              
            }
          }
          
          if (w <= Distances[v]) {
            for (int i=m_IndGr[v]; i< m_IndGr[v+1]; i++){
              
              
              int v2 = m_NodeGr[i];                                                     
              double w2 = m_WGr[i];
              
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
        
        m_result(k,i)=Dist[i];
      }
     
      
      
    }
    
  }
  
};





// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::NumericMatrix par_Bidir_mat2(std::vector<int> &dep, std::vector<int> &arr,std::vector<int> &gfrom,std::vector<int> &gto,std::vector<double> &gw,int NbNodes,std::vector<int> &Rank){
  
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
  
  int count=0;
  for (unsigned int i = 0; i < NbEdges; ++i) {
    
    if (Rank[gfrom[i]] < Rank[gto[i]])  G[gfrom[i]].push_back(std::make_pair(gto[i], gw[i]));
    if (Rank[gfrom[i]] > Rank[gto[i]]) {
      Gr[gto[i]].push_back(std::make_pair(gfrom[i], gw[i]));
      count+=1;
    }
    
    
    
  }
  
  
  //Forward
  
  std::vector<double> Distances(NbNodes, std::numeric_limits<double>::max()); 
  std::vector <int> Visited(NbNodes,0);
  
  
  for (unsigned int k=0; k!=dep.size();k++){
    if (k % 256){
      Rcpp::checkUserInterrupt ();
    }
    
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
  
  int count2=0;
  for (int i=0; i < BucketF.size();i++){
    count2+= BucketF[i].size();
    std::sort(BucketF[i].begin(),BucketF[i].end());
  }
  
  std::vector<std::vector<std::pair<int, double> > > ().swap(G);
  
  //Backward
  //Graph vectors
  Rcpp::IntegerVector NodeGr(count);
  Rcpp::NumericVector WGr(count);
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
  //Buckets
  Rcpp::IntegerVector NodeBu(count2);
  Rcpp::NumericVector WBu(count2);
  Rcpp::IntegerVector IndBu(NbNodes+1);
  count=0;
  for (int i=0; i < BucketF.size();i++){
    IndBu[i]=count;
    
    for (int j=0; j < BucketF[i].size();j++){
      NodeBu[count]=BucketF[i][j].first;
      WBu[count]=BucketF[i][j].second;
      count+=1;
    }
  }
  IndBu[NbNodes]=count;
  
  std::vector<std::vector<std::pair<int, double> > > ().swap(BucketF);
  
  Rcpp::NumericMatrix result(arr.size(),dep.size());
  

  ParBidirMat Dijfunc(NodeGr,WGr,IndGr,NodeBu,WBu,IndBu,Rcpp::wrap(arr),dep.size(),NbNodes,result);
  parallelFor(0,arr.size(),Dijfunc);
  
  
  return Rcpp::transpose(result);
  
  
}


