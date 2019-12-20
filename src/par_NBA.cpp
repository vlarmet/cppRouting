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


struct ParNBA : public Worker
{
  //input
  const RcppParallel::RVector<int> m_NodeG;
  const RcppParallel::RVector<double> m_WG;
  const RcppParallel::RVector<int> m_IndG;
  const RcppParallel::RVector<int> m_NodeGr;
  const RcppParallel::RVector<double> m_WGr;
  const RcppParallel::RVector<int> m_IndGr;
  const RcppParallel::RVector<double> m_lat;
  const RcppParallel::RVector<double> m_lon;
  const double m_k;
  RVector<int> m_dep;
  RVector<int> m_arr;
  const int m_nbnodes;
  
  //output
  RcppParallel::RVector<double> m_result;
  
  //constructor
  ParNBA(              const Rcpp::IntegerVector NodeG,
                         const Rcpp::NumericVector WG,
                         const Rcpp::IntegerVector IndG,
                         const Rcpp::IntegerVector NodeGr,
                         const Rcpp::NumericVector WGr,
                         const Rcpp::IntegerVector IndGr,
                         const Rcpp::NumericVector lat,
                         const Rcpp::NumericVector lon,
                         const double k,
                         Rcpp::IntegerVector dep,
                         Rcpp::IntegerVector arr,
                         const int nbnodes,
                         Rcpp::NumericVector result) : m_NodeG(NodeG),m_WG(WG),m_IndG(IndG),
                         m_NodeGr(NodeGr),m_WGr(WGr),m_IndGr(IndGr),m_lat(lat),m_lon(lon),m_k(k),
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
    std::vector <int> Visited1check(m_nbnodes,0);
    std::vector <int> Visited2check(m_nbnodes,0);
    
    for (std::size_t k=begin; k!=end;k++){
      
      int StartNode=m_dep[k];
      int EndNode=m_arr[k];
      
      if (StartNode==EndNode){
        m_result[k]=0;
        continue;
      }
      
      double lata=m_lat[EndNode];
      double lona=m_lon[EndNode];
      double lata2=m_lat[StartNode];
      double lona2=m_lon[StartNode];
      
      Distances[StartNode] = 0.0;  
      Visited1check[StartNode]=1;
      Distances2[EndNode] = 0.0;
      Visited2check[EndNode]=1;
      std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double> >, comp > Q;
      std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double> >, comp > Qr;
      Q.push(std::make_pair(StartNode, sqrt(pow(m_lat[StartNode]-lata,2)+pow(m_lon[StartNode]-lona,2))/m_k)); 
      Qr.push(std::make_pair(EndNode, sqrt(pow(m_lat[EndNode]-lata2,2)+pow(m_lon[EndNode]-lona2,2))/m_k)); 
      
      
      
      //double Pr=0.5*sqrt(pow(lat[StartNode]-lata,2)+pow(lon[StartNode]-lona,2))/k;
      double total1=sqrt(pow(lata2-lata,2)+pow(lona2-lona,2))/m_k;
      double total2=sqrt(pow(lata2-lata,2)+pow(lona2-lona,2))/m_k;
      int mid;
      double mu=std::numeric_limits<double>::max();
      
      while (!Q.empty() && !Qr.empty()) {  
        //Forward
        if (Q.size() < Qr.size()){
          int v=Q.top().first;
          Q.pop();
          if (Visited[v]==0){
            Visited[v]=1;
            
            if ((Distances[v] + sqrt(pow(m_lat[v]-lata,2)+pow(m_lon[v]-lona,2))/m_k) >= mu || (Distances[v] + total2 - sqrt(pow(m_lat[v]-lata2,2)+pow(m_lon[v]-lona2,2))/m_k) >= mu){}
            
            else {
              for (int i=m_IndG[v]; i< m_IndG[v+1]; i++){
                int v2 = m_NodeG[i];                                                     
                double w2 = m_WG[i];
                if (Visited[v2]==1){
                  continue;
                }
                double tentative=Distances[v]+w2;
                
                if (Visited1check[v2]==0  || Distances[v2] > tentative){
                  Distances[v2]=tentative;
                  Visited1check[v2]=1;
                  
                  Q.push(std::make_pair(v2, tentative + sqrt(pow(m_lat[v2]-lata,2)+pow(m_lon[v2]-lona,2))/m_k));
                  
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
            
            if ((Distances2[vv] + sqrt(pow(m_lat[vv]-lata2,2)+pow(m_lon[vv]-lona2,2))/m_k) >= mu || (Distances2[vv] + total1 - sqrt(pow(m_lat[vv]-lata,2)+pow(m_lon[vv]-lona,2))/m_k) >= mu){}
            
            else {
              for (int i=m_IndGr[vv]; i< m_IndGr[vv+1]; i++){
                int vv2 = m_NodeGr[i];                                                     
                double ww2 = m_WGr[i];
                if (Visited[vv2]==1){
                  continue;
                }
                double tentative=Distances2[vv]+ww2;
                
                if (Visited2check[vv2]==0  || Distances2[vv2] > tentative){
                  Distances2[vv2]=tentative;
                  Visited2check[vv2]=1;
                  
                  Qr.push(std::make_pair(vv2, tentative + sqrt(pow(m_lat[vv2]-lata2,2)+pow(m_lon[vv2]-lona2,2))/m_k));
                  
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
      
      
      if (mu>=std::numeric_limits<double>::max()){
        m_result[k] = Rcpp::NumericVector::get_na();
      }
      else {
        m_result[k]=mu;
      }
      //Reinitialize
      std::fill(Distances.begin(),Distances.end(),std::numeric_limits<double>::max());
      std::fill(Distances2.begin(),Distances2.end(),std::numeric_limits<double>::max());
      std::fill(Visited.begin(),Visited.end(),0);
      std::fill(Visited1check.begin(),Visited1check.end(),0);
      std::fill(Visited2check.begin(),Visited2check.end(),0);
      
    }
    
  }
  
};


// [[Rcpp::export]]
Rcpp::NumericVector NBA_par(Rcpp::IntegerVector dep, Rcpp::IntegerVector arr,Rcpp::IntegerVector gfrom,Rcpp::IntegerVector gto,Rcpp::NumericVector gw,int NbNodes,
                              Rcpp::NumericVector lat,Rcpp::NumericVector lon, double k){
  
  
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
  //Reverse
  
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
  
  
  
  ParNBA Dijfunc(NodeG,WG,IndG,NodeGr,WGr,IndGr,lat,lon,k,dep,arr,NbNodes,result);
  parallelFor(0,dep.length(),Dijfunc);
  
  return Rcpp::wrap(result);
  
  
  
}


