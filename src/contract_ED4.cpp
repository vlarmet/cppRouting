#include <iostream>
#include <queue>
#include <vector>
#include <fstream>
#include <limits>
#include <functional>
#include <chrono>
#include <thread>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <Rcpp.h>
#include "functions_one_to_all_ED4.h"



using namespace std;

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::List Contract_ED(std::vector<int> &gfrom,std::vector<int> &gto,std::vector<double> &gw,int NbNodes,bool display_progress){
  //Comparator 
  struct comp{
    
    bool operator()(const std::pair<int, double> &a, const std::pair<int, double> &b){
      return a.second > b.second;
    }
  };
  
  int NbEdges=gfrom.size();
  
  //Graph
  std::vector<std::vector<std::pair<int, double> > > G(NbNodes);   
  
  
  for (unsigned int i = 0; i < NbEdges; ++i) {
    
    G[gfrom[i]].push_back(std::make_pair(gto[i], gw[i]));
    
  }
  
  std::vector<double> Distances(NbNodes, std::numeric_limits<double>::max());   
  std::vector<std::vector<std::pair<int, double> > > OriginalGraph=G; 
  if (display_progress) Rcpp::Rcout << "Cleaning graph..."<<std::endl;
  for (int i=0; i< NbNodes;i++){
    std::vector<int> arr;
    std::vector<double> lim;
    for (int j=0; j < G[i].size();j++){
      arr.push_back(G[i][j].first);
      lim.push_back(G[i][j].second);
      
    }
    if (arr.size()==0) continue;
    
    Dijkstra_mod2(G,
                  OriginalGraph,
                  i,
                  arr,
                  lim,
                  NbNodes,
                  Distances);
  }
  
  //Free space
  std::vector<std::vector<std::pair<int, double> > > ().swap(OriginalGraph);
  
  //Reversed graph
  std::vector<std::vector<std::pair<int, double> > > Gr(NbNodes);
  
  
  for (unsigned int i = 0; i < G.size(); ++i) {
    for (unsigned int j = 0; j < G[i].size(); ++j) {
      Gr[G[i][j].first].push_back(std::make_pair(i,G[i][j].second));
    }
    
    
  }
  
  
  //ordering
  if (display_progress) Rcpp::Rcout << "Initial node ordering..."<< std::endl;
  priority_queue<std::pair<int, int>, vector<std::pair<int, int> >, comp > Queue;
  
  std::vector<int> Contracted(NbNodes,0);
  
  
  
  for (int i=0; i < NbNodes; i++){
    Rcpp::checkUserInterrupt();
    
    
    int ED=Edge_dif(i,G,Gr,NbNodes,Distances);
    
    
    Queue.push(std::make_pair(i,ED));
  }
  
  //Contracting
  //Final augmented graphs
  std::vector<std::vector<std::pair<int, double> > > AugG=G;   
  std::vector<std::vector<std::pair<int, double> > > AugGr=Gr;
  
  //Shortcuts will be stored in 3 vectors (from, to, contracted node)
  std::vector<int> ShortF;
  std::vector<int> ShortT;
  std::vector<int> ShortC;
  
  //Free space
  std::vector<std::vector<std::pair<int, double> > >().swap(Gr);
  
  //Node order
  std::vector<int> Order(NbNodes,0);
  
  int count=0;
  int count2=0;
  int count3=0;
  bool error=false;
  if (display_progress) Rcpp::Rcout << "Contracting nodes..."<< std::endl;
  Progress p(NbNodes, display_progress);
  
  //Contracting all nodes
  while(!Queue.empty()){
    
    if (error) break;
    
    Rcpp::checkUserInterrupt();
    auto t1 = std::chrono::high_resolution_clock::now();
    
    int v=Queue.top().first;
    int imp=Queue.top().second;
    Queue.pop();
    
    if (Contracted[v]==1) {
      count3+=1;
      continue; //already contracted
    }
    
    
    int DN=0;
    for (int i=0; i < G[v].size();i++){
      if (Contracted[G[v][i].first]==1) DN+=1;
    }

    
    int New_imp=Edge_dif(v,AugG,AugGr,NbNodes,Distances)*190+DN*120;
    
    
    if (Queue.empty()){
      count+=1;
      p.increment();
      
      Order[v]=count;
      Contracted[v]=1;
      auto t2 = std::chrono::high_resolution_clock::now();
      
      Contract(v,
               AugG,
               AugGr,
               G,
               NbNodes,
               Distances,
               Contracted,count,
               error,
               ShortF,
               ShortT,
               ShortC);
    }
    
    else{
      if (New_imp > Queue.top().second){ //compare to the next min
        count2+=1;
        //reinsert
        Queue.push(std::make_pair(v,New_imp));
        
        
        
        
      }
      
      
      else{
        
        //contract
        
        count+=1;
        p.increment();
        Order[v]=count;
        Contracted[v]=1;
        
        Contract(v,
                 AugG,
                 AugGr,
                 G,
                 NbNodes,
                 Distances,
                 Contracted,count,
                 error,
                 ShortF,
                 ShortT,
                 ShortC);
        
        
      }
      
    }
    
  }
  
  
  //free space
  std::vector<std::vector<std::pair<int, double> > >().swap(AugG);
  std::vector<std::vector<std::pair<int, double> > >().swap(AugGr);
  
  
  //Augmented graph
  count=0;
  for (int i=0; i < NbNodes;i++) count+=G[i].size();
  
  
  std::vector<int> Newfrom(count);
  std::vector<int> Newto(count);
  std::vector<double> Neww(count);
  
  int index=0;
  for (int i=0; i < NbNodes;i++){
    for (int j=0; j < G[i].size();j++){
      
      Newfrom[index]=i;
      Newto[index]= G[i][j].first;
      Neww[index]= G[i][j].second;
      index+=1;
      
    }
  }
  
  //Return
  
  Rcpp::List final(3);
  final[0]=Newfrom;
  final[1]=Newto;
  final[2]=Neww;
  
  Rcpp::List Shortcuts(3);
  Shortcuts[0]=ShortF;
  Shortcuts[1]=ShortT;
  Shortcuts[2]=ShortC;
  
  //Rcpp::NumericVector result(NbNodes);
  Rcpp::List finalList(3);
  finalList[0]=final;
  finalList[1]=Order;
  finalList[2]=Shortcuts;
  return (finalList); 
}