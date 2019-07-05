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

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::List Simplify2(std::vector<int> gfrom,std::vector<int> gto,std::vector<double> gw,int NbNodes,bool loop,std::vector<int> keep,std::vector<std::string> dict){
  
  //Rcpp::IntegerVector Remove(NbNodes,0);
  
  //Graphs
  
  int NbEdges=gfrom.size();
  
  std::vector<std::vector<int> > G(NbNodes);   
  std::vector<std::vector<int> > Gr(NbNodes);
  std::vector<std::vector<std::pair<int,double> > > Graph(NbNodes);   
  //std::vector<std::vector<std::pair<int,double> > > Graphr(NbNodes);
  
  for (unsigned int i = 0; i < NbEdges; ++i) {
    
    G[gfrom[i]].push_back(gto[i]);
    Gr[gto[i]].push_back(gfrom[i]);
    Graph[gfrom[i]].push_back(std::make_pair(gto[i],gw[i]));
  }
  

  std::vector<int> Junction(NbNodes,0);
  
  // Detect junction nodes
  for (int i=0; i < G.size() ; ++i){
    //Nodes to keep
    if (keep[i]==1){
      continue;
    }
    
    //undirected
    if ((G[i].size()==2 && Gr[i].size()==2)){
      std::vector<int> vec1=G[i];
      std::vector<int> vec2=Gr[i];
      std::sort(vec1.begin(),vec1.end());
      std::sort(vec2.begin(),vec2.end());
      
      if (vec1[0]==vec2[0] && vec1[1]==vec2[1]){
        Junction[i]=2;
       
      }
      
      
    }
    
    //directed
    if (G[i].size()==1 && Gr[i].size()==1){
      if (G[i]!=Gr[i]){
        Junction[i]=1;
        
        
      } 
      
      
    }
    
  }
  
  //Edge index
  std::vector<std::vector<int> > Locfrom(NbNodes); //node, edge index 
  std::vector<std::vector<int> > Locto(NbNodes);
  
  for (int i=0; i < gfrom.size(); ++i){
    if (Junction[gfrom[i]]>0){
      Locfrom[gfrom[i]].push_back(i);
    }
    if (Junction[gto[i]]>0){
      Locto[gto[i]].push_back(i);
    }
  }
  
  /// Construct undirected edges 
  
  
  
  std::vector<std::vector<int> > Edges;
  std::vector<int> Treated(NbNodes,0);
  
  for (int i=0; i < NbNodes ; ++i){
    if (Junction[i]!=2 || Treated[i]==1){
      continue;
    }
    
    
    
    //Rcpp::Rcerr << i << std::endl;
    Rcpp::checkUserInterrupt();
    
    std::vector<int> path;
    int node1=i;
    Treated[i]=1;
    path.push_back(i);
    
    

    
   if (Locto[node1].size()!=2){
     continue;
   }
   
   int breaked=0;
   
   for (int j=0; j<2; ++j){
     int node2=gfrom[Locto[node1][j]];
     
     breaked=0;
     path.insert(path.begin(),node2);
     Treated[node2]=1;
     
     while (Junction[node2]==2){
       Rcpp::checkUserInterrupt();
       
       
       int node=gfrom[Locto[node2][0]];
       
       
       if (Treated[node]==1 && Junction[node]==2){

          node=gfrom[Locto[node2][1]];
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
    for (int i=0; i < NbNodes ; ++i){
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
      
      if (Locto[node1].size()!=1){
        continue;
      }
      
      int node2=gfrom[Locto[node1][0]];
      
      path.insert(path.begin(),node2);
      Treated[node2]=1;
      
      //One direction
      while (Junction[node2]==1){
        Rcpp::checkUserInterrupt();
        
        int node=gfrom[Locto[node2][0]];

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
      
      node2=gto[Locfrom[node1][0]];
      path.push_back(node2);
      Treated[node2]=1;
      while (Junction[node2]==1){
        Rcpp::checkUserInterrupt();
        
        int node=gto[Locfrom[node2][0]];
        
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
    
    
  
  

  //Rcpp::Rcerr <<"ok"<<std::endl;
  // Reconstruct symbolic graph with contracted edges
  std::vector<int> Remove(NbEdges,0); //edge index to remove
  std::vector<int> Removenode(NbNodes,0);  // nodes inserted in new edge
  
    
  std::vector<double> Acc;
  std::vector<double> Acc2;
  for (int i=0; i < Edges.size(); ++i){
    std::vector<int> path=Edges[i];
    
    for (int k=1; k < (path.size()-1); ++k){
      
      for (int it=0; it < Locto[path[k]].size(); ++it){
        Remove[Locto[path[k]][it]]=1;
        Remove[Locfrom[path[k]][it]]=1;
        
        Removenode[path[k]]=1;
      }
      
    }
    
    double w=0.0;
    
    for (int k=0; k < (path.size()-1); ++k){
      int node1=path[k];
      int node2=path[k+1];
      
      for (int m=0; m < Graph[node1].size(); ++m){
        if (Graph[node1][m].first==node2){
          w+=Graph[node1][m].second;
        }
      }
      
      
      
    }
    Acc.push_back(w);
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
    Acc2.push_back(w2);

  }
  
  //Remove loops
  if (loop==true){
    for (int i=0; i < NbNodes; ++i){
      if (Junction[i]>0 && Removenode[i]==0){
        for (int it=0; it < Locto[i].size(); ++it){
          Remove[Locto[i][it]]=1;
          Remove[Locfrom[i][it]]=1;
        }
      }
      
    }
    
  }
  
  
  
  //Rcpp::Rcerr <<"ok"<<std::endl;

    std::vector<int> Newfrom;
    std::vector<int> Newto;
    std::vector<double> Neww;
    
    for (int i=0; i < gfrom.size(); ++i){
      if (Remove[i]==0){
        Newfrom.push_back(gfrom[i]);
        Newto.push_back(gto[i]);
        Neww.push_back(gw[i]);
      }
      
      
    }
    //Rcpp::Rcerr <<"ok"<<std::endl;
    
    std::vector<std::vector<int> > validEdge;//valid edges to keep (no loops)
    for (int i=0; i < Edges.size(); ++i){
      //avoid construct new loops
      if (Edges[i][0]==Edges[i][Edges[i].size()-1] && keep[Edges[i][0]]==0){
        continue;
      }
      validEdge.push_back(Edges[i]);

      if (Junction[Edges[i][1]]==2){
        Newfrom.push_back(Edges[i][0]);
        Newto.push_back(Edges[i][Edges[i].size()-1]);
        Neww.push_back(Acc[i]);
        Newfrom.push_back(Edges[i][Edges[i].size()-1]);
        Newto.push_back(Edges[i][0]);
        Neww.push_back(Acc2[i]);
        
      }
      if (Junction[Edges[i][1]]==1) {
        Newfrom.push_back(Edges[i][0]);
        Newto.push_back(Edges[i][Edges[i].size()-1]);
        Neww.push_back(Acc[i]);
        
      }
      
    } 
    

    
    // Number of unique nodes
    std::unordered_set<int> nodes;
    for (int i=0; i < Newfrom.size(); ++i){
      nodes.insert(Newfrom[i]);
      nodes.insert(Newto[i]);
    }
    
    //Convert node id to string in edge list
    
    std::vector<std::vector<std::string > > edge_dict(validEdge.size());
    for (int i=0; i < validEdge.size(); ++i){
      std::vector<std::string> temp(validEdge[i].size());
      for (int j=0; j < validEdge[i].size(); ++j){
        temp[j]=dict[validEdge[i][j]];
      }
      
      edge_dict[i]=temp;
    }
    
  

  Rcpp::NumericMatrix final(Newfrom.size(),3);
  
  
  Rcpp::IntegerVector temp=wrap(Newfrom);
  final( _ ,0)=temp;
  Rcpp::IntegerVector temp2=wrap(Newto);
  final( _ ,1)=temp2;
  Rcpp::NumericVector temp3=wrap(Neww);
  final( _ ,2)=temp3;
  
  //return Rcpp::wrap(Edges);
  
  Rcpp::List finalList(3);
  finalList[0]=final;
  finalList[1]=nodes.size();
  Rcpp::List edgelist=wrap(edge_dict);
  finalList[2]=edgelist;
  
  
  
  return finalList;
   
    
    
  
  
  
  
  
}