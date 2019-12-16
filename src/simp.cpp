#include <iostream>
#include <queue>
#include <vector>
#include <fstream>
#include <limits>
#include <functional>
#include <algorithm>
#include <string>
#include <Rcpp.h>




using namespace Rcpp;

void quickDelete( int idx ,std::vector<std::pair<int,double> > &vec)
{
  vec[idx] = vec.back();
  vec.pop_back();
}

void quickDelete2( int idx ,std::vector<int> &vec)
{
  vec[idx] = vec.back();
  vec.pop_back();
}


void Simplify(std::vector<std::vector<std::pair<int,double> > > &Graph,
              std::vector<std::vector<int> > &Gr,
              std::vector<int> &Junction,
              std::vector<std::vector<int> > &Edges,
              std::vector<int> &keep,
              std::vector<int> &Treated,
              bool &loop,
              int &N){
  
  // Detect junction nodes
  for (int i=0; i < Graph.size() ; ++i){
    //Nodes to keep
    if (keep[i]==1){
      continue;
    }
    
    //undirected
    if ((Graph[i].size()==2 && Gr[i].size()==2)){
      std::vector<std::pair<int,double> > vec1=Graph[i];
      std::vector<int> vec2=Gr[i];
      std::sort(vec1.begin(),vec1.end());
      std::sort(vec2.begin(),vec2.end());
      
      if (vec1[0].first==vec2[0] && vec1[1].first==vec2[1]){
        Junction[i]=2;
        
      }
      
      
    }
    
    //directed
    if (Graph[i].size()==1 && Gr[i].size()==1){
      if (Graph[i][0].first!=Gr[i][0]){
        Junction[i]=1;
        
        
      } 
      
      
    }
    
  }
  ////////////////////////////////////////////////////////////////////////////////
  for (int i=0; i < Graph.size() ; ++i){
    if (Junction[i]!=2 || Treated[i]==1){
      continue;
    }
    
    
    
    //Rcpp::Rcerr << i << std::endl;
    Rcpp::checkUserInterrupt();
    
    std::vector<int> path;
    int node1=i;
    Treated[i]=1;
    path.push_back(i);
    
    
    
    
    if (Gr[node1].size()!=2){
      continue;
    }
    
    int breaked=0;
    
    for (int j=0; j<2; ++j){
      int node2=Gr[node1][j];
      
      breaked=0;
      path.insert(path.begin(),node2);
      Treated[node2]=1;
      
      while (Junction[node2]==2){
        Rcpp::checkUserInterrupt();
        
        
        int node=Gr[node2][0];
        
        
        if (Treated[node]==1 && Junction[node]==2){
          
          node=Gr[node2][1];
          //loop in the graph
          if (Treated[node]==1 && Junction[node]==2){
            breaked=1;
            break;
          }
          
        }
        //Rcpp::Rcerr << node << std::endl;
        
        node2=node;
        
        //loop in the graph
        
        if (node==i){
          breaked=1;
          break;
        }
        
        path.insert(path.begin(),node2);
        Treated[node2]=1;
        
      }
      
      std::reverse(path.begin(),path.end());
      
      
    }
    
    
    
    // for (int k=0; k < path.size();++k){
    //  Rcpp::Rcerr << path[k] << "__";
    //}
    
    
    if (breaked==0){
      
      
      Edges.push_back(path);
    }
    
    
    
    
  }
  
  
  int breaked=0;
  /// Construct directed edges 
  for (int i=0; i < Graph.size() ; ++i){
    if (Junction[i]!=1 || Treated[i]==1){
      continue;
    }
    
    breaked=0;
    
    
    //Rcpp::Rcerr << i << std::endl;
    Rcpp::checkUserInterrupt();
    
    std::vector<int> path;
    int node1=i;
    path.push_back(i);
    Treated[i]=1;
    
    if (Gr[node1].size()!=1){
      continue;
    }
    
    int node2=Gr[node1][0];
    
    path.insert(path.begin(),node2);
    Treated[node2]=1;
    
    //One direction
    while (Junction[node2]==1){
      Rcpp::checkUserInterrupt();
      
      int node=Gr[node2][0];
      
      //Rcpp::Rcerr << node << std::endl;
      
      node2=node;
      //loop in the graph
      if (node==i){
        breaked=1;
        break;
      }
      
      
      path.insert(path.begin(),node2);
      Treated[node2]=1;
      
    }
    
    //other direction
    
    node2=Graph[node1][0].first;
    path.push_back(node2);
    Treated[node2]=1;
    while (Junction[node2]==1){
      Rcpp::checkUserInterrupt();
      
      int node=Graph[node2][0].first;
      
      //Rcpp::Rcerr << node << std::endl;
      
      node2=node;
      //loop in the graph
      if (node==i){
        breaked=1;
        break;
      }
      
      
      path.push_back(node2);
      Treated[node2]=1;
      
    }
    
    
    if (breaked==0){
      Edges.push_back(path);
    }
    
    
    
    
  }
  
  /////////////////////////////////////////////////////////////////////////////////
  std::vector<int> Removenode(Graph.size(),0);  // nodes inserted in new edge
  
  
  std::vector<double> Acc(Edges.size());
  std::vector<double> Acc2(Edges.size());
  for (int i=0; i < Edges.size(); ++i){
    std::vector<int> path=Edges[i];
    double w=0.0;
    for (int k=0; k < (path.size()-1); ++k){
      
      Removenode[path[k]]=1;
      int node1=path[k];
      int node2=path[k+1];
      
      for (int m=0; m < Graph[node1].size(); ++m){
        if (Graph[node1][m].first==node2){
          w+=Graph[node1][m].second;
        }
      }
      
    }
    
    
    

    Acc[i]=w;
    //Autre sens
    std::vector<int> path2=path;
    std::reverse(path2.begin(),path2.end());
    double w2=0.0;
    for (int k=0; k < (path2.size()-1); ++k){
      int node1=path2[k];
      int node2=path2[k+1];
      
      for (int m=0; m < Graph[node1].size(); ++m){
        if (Graph[node1][m].first==node2){
          w2+=Graph[node1][m].second;
        }
      }
      
      
      
    }
    Acc2[i]=w2;
    
  }
  
  ///////////////////////////////////////////////////////////////////////
  //Remove from graph and gr
  for (int i=0; i < Edges.size();i++){
    for (int j=1; j < Edges[i].size()-1;j++){
      //Graph
      std::vector<std::pair<int,double > > Gr_remov=Graph[Edges[i][j]];
      
      Graph[Edges[i][j]].erase(Graph[Edges[i][j]].begin(),Graph[Edges[i][j]].end());
      for (int k=0; k < Gr[Edges[i][j]].size();k++){
        for (int m=0; m < Graph[Gr[Edges[i][j]][k]].size();m++){
          if (Graph[Gr[Edges[i][j]][k]][m].first==Edges[i][j]) {
            quickDelete(m,Graph[Gr[Edges[i][j]][k]]);
          }
        }
        
      }
      
      //Gr
      Gr[Edges[i][j]].erase(Gr[Edges[i][j]].begin(),Gr[Edges[i][j]].end());
      for (int k=0; k < Gr_remov.size();k++){
        for (int m=0; m < Gr[Gr_remov[k].first].size();m++){
          if (Gr[Gr_remov[k].first][m]==Edges[i][j]){
            
            quickDelete2(m,Gr[Gr_remov[k].first]);
          }
        }
      }
    }
  }
  
  
  //Remove loops
  if (loop==true){
    for (int i=0; i < Graph.size(); ++i){
      if (Junction[i]>0 && Removenode[i]==0){
        std::vector<std::pair<int,double > > Gr_remov=Graph[i];
        
        Graph[i].erase(Graph[i].begin(),Graph[i].end());
        for (int k=0; k < Gr[i].size();k++){
          for (int m=0; m < Graph[Gr[i][k]].size();m++){
            if (Graph[Gr[i][k]][m].first==i) {
              quickDelete(m,Graph[Gr[i][k]]);
            }
          }
          
        }
        
        Gr[i].erase(Gr[i].begin(),Gr[i].end());
        for (int k=0; k < Gr_remov.size();k++){
          for (int m=0; m < Gr[Gr_remov[k].first].size();m++){
            if (Gr[Gr_remov[k].first][m]==i){
              
              quickDelete2(m,Gr[Gr_remov[k].first]);
            }
          }
        }
      }
      
    }
    
  }
  
  //Add shortcuts and avoid duplicated edges
  bool dup;
  for (int i=0; i < Edges.size(); i++){
    if (Edges[i][0]==Edges[i][Edges[i].size()-1] && keep[Edges[i][0]]==0){
      continue;
    }
    
    
    dup=false;
    int node1=Edges[i][0];
    int node2=Edges[i][Edges[i].size()-1];
    
    for (int k=0; k<Graph[node1].size(); k++){
      if (Graph[node1][k].first==node2){
        dup=true;
        if (Graph[node1][k].second >= Acc[i]){
          Graph[node1][k].second=Acc[i];
          break;
        }
        
        
      }
    }
    
    if (dup==false) {
      Graph[node1].push_back(std::make_pair(node2,Acc[i]));
      Gr[node2].push_back(node1);
    }
    
    //Other direction
    if (Junction[Edges[i][1]]==2){
      dup=false;
      for (int k=0; k<Graph[node2].size(); k++){
        if (Graph[node2][k].first==node1){
          dup=true;
          if (Graph[node2][k].second >= Acc2[i]){
            Graph[node2][k].second=Acc2[i];
            break;
          }
          
          
        }
      }
      
      if (dup==false) {
        Graph[node2].push_back(std::make_pair(node1,Acc2[i]));
        Gr[node1].push_back(node2);
      }
    }
    
    
  }
  
  
  ///////////////////////////////////////////////////////////////////
  //Number of nodes
  std::vector<int> Nodes(Graph.size(),0);
  for (int i=0; i < Graph.size();i++){
    for (int j=0; j < Graph[i].size();j++){
      Nodes[Graph[i][j].first]+=1;
    }
    
    for (int j=0; j < Gr[i].size();j++){
      Nodes[Gr[i][j]]+=1;
    }
    
  }
  
  int count=0;
  for (int i=0; i < Nodes.size();i++){
    if (Nodes[i]>0) count+=1;
  }
  
  N=count;
  
  ////////////////////////////////////////////////////////////////
  //Reinitialize
  std::fill(Treated.begin(),Treated.end(),0);
  std::fill(Junction.begin(),Junction.end(),0);
  Edges.resize(0);
  
  
  
  
  
}