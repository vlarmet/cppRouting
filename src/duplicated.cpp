#include <iostream>
#include <queue>
#include <vector>
#include <fstream>
#include <limits>
#include <functional>
#include <algorithm>
#include <Rcpp.h>




using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::DataFrame Remove_duplicate(std::vector<int> gfrom,std::vector<int> gto,std::vector<double> gw,int NbNodes){

  std::vector<std::vector<std::pair<int,double> > > Graph(NbNodes);   
  //std::vector<std::vector<std::pair<int,double> > > Graphr(NbNodes);
  int count=0;
  for (unsigned int i = 0; i < gfrom.size(); ++i) {

    count=0;
    if (Graph[gfrom[i]].size()==0){
      Graph[gfrom[i]].push_back(std::make_pair(gto[i],gw[i]));
    }
    else {
    
    for (int j=0;j < Graph[gfrom[i]].size();++j){
      
      if (Graph[gfrom[i]][j].first==gto[i]){
        if (Graph[gfrom[i]][j].second >= gw[i]){
          count=1;
          Graph[gfrom[i]][j].second=gw[i];
        }
        else {
          count=1;
          break;
        }
    
      }

    }
    
    if (count==0){
      Graph[gfrom[i]].push_back(std::make_pair(gto[i],gw[i]));
    }
    }
    
  }
  
  std::vector<int> From;
  std::vector<int> To;
  std::vector<double> W;
  
  for (int i=0; i < Graph.size(); ++i){
    for (int j=0; j < Graph[i].size(); ++j){
      From.push_back(i);
      To.push_back(Graph[i][j].first);
      W.push_back(Graph[i][j].second);
    }
    
  }
  
  Rcpp::DataFrame final=Rcpp::DataFrame::create(Rcpp::Named("from")=From,
                                                Rcpp::Named("to")=To,
                                                Rcpp::Named("dist")=W);

  
  return final;
  
  
  
}
