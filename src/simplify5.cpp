#include <iostream>
#include <queue>
#include <vector>
#include <fstream>
#include <limits>
#include <functional>
#include <algorithm>
#include <string>
#include "simp.h"
#include <Rcpp.h>




using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::List Simplify3(std::vector<int> gfrom,std::vector<int> gto,std::vector<double> gw,
                     int NbNodes,bool loop,std::vector<int> keep,
                     bool iterate, bool progress){
  
  //Rcpp::IntegerVector Remove(NbNodes,0);
  
  //Graphs
  
  
  int NbEdges=gfrom.size();
  
  std::vector<std::vector<std::pair<int,double> > > Graph(NbNodes);   
  std::vector<std::vector<int> > Gr(NbNodes);
  //std::vector<std::vector<std::pair<int,double> > > Graphr(NbNodes);
  
  for (unsigned int i = 0; i < NbEdges; ++i) {
    
    Gr[gto[i]].push_back(gfrom[i]);
    Graph[gfrom[i]].push_back(std::make_pair(gto[i],gw[i]));
  }
  
  std::vector<int> Treated(NbNodes,0);
  std::vector<int> Junction(NbNodes,0);
  std::vector<std::vector<int> > Edges;
  
  
  int count=0;
  
  int N=NbNodes;
  
  if (iterate){
    while (true){
      
      int N2=N;
      Simplify(Graph,Gr,Junction,Edges,keep,Treated,loop,N);
      count+=1;
      if (progress && (N2-N)!=0) Rcpp::Rcout << "iteration : "<<count<< " - " <<N2-N<< " nodes removed" << std::endl;
      if (N2==N) break;
      
    }
    
  }
  else{
    int N2=N;
    Simplify(Graph,Gr,Junction,Edges,keep,Treated,loop,N);
    
    if (progress) Rcpp::Rcout<< N2-N << " nodes removed" << std::endl;
  }
  
  
  
  
  
  std::vector<int> from2;
  std::vector<int> to2;
  std::vector<double> w2;
  for (int i=0; i < Graph.size();i++){
    for (int j=0; j < Graph[i].size();j++){
      from2.push_back(i);
      to2.push_back(Graph[i][j].first);
      w2.push_back(Graph[i][j].second);
    }
    
  }

  Rcpp::List finalList(3);
  finalList[0]=from2;
  finalList[1]=to2;
  finalList[2]=w2;

  
  return finalList;
  
  
  
  
  
  
  
  
}