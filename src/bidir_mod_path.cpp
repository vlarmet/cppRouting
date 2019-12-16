#include <iostream>
#include <queue>
#include <vector>
#include <fstream>
#include <limits>
#include <functional>
#include <chrono>
#include <thread>
#include <string>
#include <Rcpp.h>




using namespace std;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::List Bidir_mod_path(std::vector<int> &dep, std::vector<int> &arr,std::vector<int> &gfrom,std::vector<int> &gto,std::vector<double> &gw,int NbNodes,
                          std::vector<int> &Rank,std::vector<std::string> dict,
                          std::vector<int> &keep,
                          std::vector<int> ShortF,std::vector<int> ShortT,std::vector<int> ShortC){
  
  //Shortcuts
  std::vector<std::vector<std::pair<int,int> > > Shortcuts(NbNodes);
  
  for (int i=0; i < ShortF.size();i++){
    Shortcuts[ShortF[i]].push_back(std::make_pair(ShortT[i],ShortC[i]));
    
  }
  
  std::vector<std::vector<std::string> > result(dep.size());
  
  
  
  
  struct comp{
    
    bool operator()(const std::pair<int, double> &a, const std::pair<int, double> &b){
      return a.second > b.second;
    }
  };
  
  //Graphs
  
  int NbEdges=gfrom.size();
  
  std::vector<std::vector<std::pair<int, double> > > G(NbNodes);   
  std::vector<std::vector<std::pair<int, double> > > Gr(NbNodes);
  
  for (unsigned int i = 0; i < NbEdges; ++i) {
    
    
    if (Rank[gfrom[i]] < Rank[gto[i]])  G[gfrom[i]].push_back(std::make_pair(gto[i], gw[i]));
    if (Rank[gfrom[i]] > Rank[gto[i]]) Gr[gto[i]].push_back(std::make_pair(gfrom[i], gw[i]));
    
    
  }
  
  
  //Boucle sur chaque trajet
  
  std::vector<double> Distances(NbNodes, std::numeric_limits<double>::max()); 
  std::vector<double> Distances2(NbNodes, std::numeric_limits<double>::max()); 
  std::vector <int> Visited1(NbNodes,0);
  std::vector <int> Visited2(NbNodes,0);
  std::vector <int> Visited;
  std::vector<int> Parents(NbNodes, -1);     
  std::vector<int> Parents2(NbNodes, -1);   
  
  
  
  for (unsigned int k=0; k!=dep.size();k++){
    if (k % 256){
      Rcpp::checkUserInterrupt ();
    }
    
    auto t1 = std::chrono::high_resolution_clock::now();
    int StartNode=dep[k];
    int EndNode=arr[k];
    
    
    
    
    Distances[StartNode] = 0.0;  
    Distances2[EndNode] = 0.0;
    
    
    
    
    
    std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double> >, comp > Q;
    std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double> >, comp > Qr;
    Q.push(std::make_pair(StartNode, 0.0)); 
    Qr.push(std::make_pair(EndNode, 0.0));
    
    int mid;
    double mu=std::numeric_limits<double>::max();
    
    

    
    
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
          mid=v;
          mu=Distances[v]+Distances2[v];
          
        }
        
        
        
        
        
        if (w <= Distances[v]) {
          for (int i=0; i< G[v].size(); i++){
            int v2 = G[v][i].first;                                                     
            double w2 = G[v][i].second;

            
            if (Distances[v] + w2 < Distances[v2]) {                               
              Distances[v2] = Distances[v] + w2;                                   
              Q.push(std::make_pair(v2, Distances[v2]));
              Visited1[v2]=1;
              Parents[v2] = v;
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
          mid=vv;
          mu=Distances[vv]+Distances2[vv];
          
        }
        
        if (ww <= Distances2[vv]) {
          for (int i=0; i< Gr[vv].size(); i++){
            int vv2 = Gr[vv][i].first;                                                     
            double ww2 = Gr[vv][i].second;
            
            if (Distances2[vv] + ww2 < Distances2[vv2]) {                               
              Distances2[vv2] = Distances2[vv] + ww2;                                   
              
              Qr.push(std::make_pair(vv2, Distances2[vv2]));
              Visited2[vv2]=1;
              Parents2[vv2] = vv;
              Visited.push_back(vv2);
              
            }
          }
        }
        
      }
      
    }
    
    
    

    
    //Extracting path
    std::vector <int> result2;
    if (mu<std::numeric_limits<double>::max() ){
      for (auto p = Parents2[mid]; p != -1; p = Parents2[p]){
        result2.insert(result2.begin(),p);
      }
      
      if (Distances[mid]!=std::numeric_limits<double>::max() || Distances2[mid]!=std::numeric_limits<double>::max()){
        result2.push_back(mid);
      }
      
      for (auto p = Parents[mid]; p != -1; p = Parents[p]){
        result2.push_back(p);
      }
      
      std::reverse(result2.begin(),result2.end());
      
      
      std::vector<int> temp(result2);
      
      
      if (result2.size() > 1){
        while (true){
          Rcpp::checkUserInterrupt();
          
          int again=0;
          std::vector<std::pair<int,int> > to_insert;
          int count=0;
          for (int i=0; i < (result2.size() -1);i++){
            
            for (int j=0;j<Shortcuts[result2[i]].size();j++){
              if (Shortcuts[result2[i]][j].first==result2[i+1]){
                to_insert.push_back(std::make_pair(i+count+1, Shortcuts[result2[i]][j].second));
                count+=1;
                again=1;
              } 
            }
            
          }
          for (int i=0; i < count; i++) temp.insert(temp.begin() + to_insert[i].first, to_insert[i].second);
          
          result2=temp;
          if (again==0) break;
          
        }
      }
    }
    
    
    for (int i=0; i<Visited.size();i++){
      Visited1[Visited[i]]=0;
      Visited2[Visited[i]]=0;
      Distances[Visited[i]]=std::numeric_limits<double>::max();
      Distances2[Visited[i]]=std::numeric_limits<double>::max();
      Parents[Visited[i]]= -1;
      Parents2[Visited[i]]= -1;
    }
    Visited.clear();
    


    std::reverse(result2.begin(),result2.end());
    std::vector<std::string > result3;
    for (int i=0; i < result2.size();i++) if (keep[result2[i]]==1) result3.push_back(dict[result2[i]]);
    

      
      
      result[k]=result3;
    
    
    
    
    
  }
  
  
  return Rcpp::wrap(result);
  
}


