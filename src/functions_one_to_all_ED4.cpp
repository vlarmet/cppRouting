#include <iostream>
#include <queue>
#include <vector>
#include <fstream>
#include <limits>
#include <functional>
#include <Rcpp.h>
#include <chrono>


using namespace std;

// [[Rcpp::plugins(cpp11)]]

void quickDelete3( int idx ,std::vector<std::pair<int,double> > &vec)
{
  vec[idx] = vec.back();
  vec.pop_back();
}



void Dijkstra_mod(std::vector<std::vector<std::pair<int, double> > > &Graph,std::vector<std::vector<std::pair<int, double> > > &Graphr,
                  std::vector<std::vector<std::pair<int, double> > > &OrGraph,std::vector<std::vector<std::pair<int, double> > > &OrGraphr,
                  int dep, std::vector<int> &arr, std::vector<double> &lim, int NbNodes, std::vector<double> &Distances){
  
  
  double maxlim=*max_element(lim.begin(), lim.end());
  
  
  
  //Comparator 
  struct comp{
    
    bool operator()(const std::pair<int, double> &a, const std::pair<int, double> &b){
      return a.second > b.second;
    }
  };
  
  
  
  //DIjkstra
  
  
  
  
  Distances[dep] = 0.0;                                                     
  
  
  priority_queue<std::pair<int, double>, vector<std::pair<int, double> >, comp > Q;
  
  
  Q.push(std::make_pair(dep, 0.0));   
  
  
  //int count=0;
  std::vector<int> visited;
  //std::vector<int>::iterator it;
  while (!Q.empty()) {                                                        
    int v = Q.top().first;                                                     
    double w = Q.top().second;                                                     
    Q.pop();
    visited.push_back(v);
    
    if (Distances[v] > maxlim){
      break;
    }
    
    
    if (w <= Distances[v]) {                                                    
      
      for (unsigned int i=0; i< OrGraph[v].size(); i++) {
        std::pair<int,double> j = OrGraph[v][i];                                               
        int v2 = j.first;                                                      
        double w2 = j.second;
        
        if (Distances[v] + w2 < Distances[v2]) {                               
          Distances[v2] = Distances[v] + w2;                                   
          
          
          Q.push(make_pair(v2, Distances[v2]));
          visited.push_back(v2);
        }
      }
    }
    
  }
  
  
  for (int i=0; i < arr.size(); i++){
    
    if (Distances[arr[i]]< lim[i]) {
      
      Graph[dep][i].second=Distances[arr[i]];
      int ind;
      
      for (int j=0; j < Graphr[arr[i]].size();j++){
        if (Graphr[arr[i]][j].first==dep) ind=j;
        break;
      }
      Graphr[arr[i]][ind].second=Distances[arr[i]];
    }
  }
  
  //Reinitialize distances vector
  for (int i=0; i < visited.size();i++) Distances[visited[i]]=std::numeric_limits<double>::max();
  
  
  
}




void Dijkstra_mod2(std::vector<std::vector<std::pair<int, double> > > &Graph,
                   std::vector<std::vector<std::pair<int, double> > > &OrGraph,
                   int dep, std::vector<int> &arr, std::vector<double> &lim, int NbNodes, std::vector<double> &Distances){
  
  
  double maxlim=*max_element(lim.begin(), lim.end());
  
  
  
  //Comparator 
  struct comp{
    
    bool operator()(const std::pair<int, double> &a, const std::pair<int, double> &b){
      return a.second > b.second;
    }
  };
  
  
  
  //DIjkstra
  
  
  
  
  Distances[dep] = 0.0;                                                     
  
  
  priority_queue<std::pair<int, double>, vector<std::pair<int, double> >, comp > Q;
  
  
  Q.push(std::make_pair(dep, 0.0));   
  
  
  //int count=0;
  std::vector<int> visited;
  //std::vector<int>::iterator it;
  while (!Q.empty()) {                                                        
    int v = Q.top().first;                                                     
    double w = Q.top().second;                                                     
    Q.pop();
    visited.push_back(v);
    
    if (Distances[v] > maxlim){
      break;
    }
    
    
    if (w <= Distances[v]) {                                                    
      
      for (unsigned int i=0; i< OrGraph[v].size(); i++) {
        std::pair<int,double> j = OrGraph[v][i];                                               
        int v2 = j.first;                                                      
        double w2 = j.second;
        
        if (Distances[v] + w2 < Distances[v2]) {                               
          Distances[v2] = Distances[v] + w2;                                   
          
          
          Q.push(make_pair(v2, Distances[v2]));
          visited.push_back(v2);
        }
      }
    }
    
  }
  
  
  for (int i=0; i < arr.size(); i++){
    
    if (Distances[arr[i]]< lim[i]) {
      
      for (int j=0; j < Graph[dep].size();j++){
        if (Graph[dep][j].first==arr[i]){
          quickDelete3(j,Graph[dep]);
          break;
        }
      }
      
      
    }
  }
  
  //Reinitialize distances vector
  for (int i=0; i < visited.size();i++) Distances[visited[i]]=std::numeric_limits<double>::max();
  
  
  
}







int Dijkstra_bool(std::vector<std::vector<std::pair<int, double> > > &Graph,int dep, std::vector<int> &arr, std::vector<double> &lim, int NbNodes, int node, std::vector<double> &Distances){
  
  
  double maxlim=*max_element(lim.begin(), lim.end());
  
  
  auto t1 = std::chrono::high_resolution_clock::now();
  //Comparator 
  struct comp{
    
    bool operator()(const std::pair<int, double> &a, const std::pair<int, double> &b){
      return a.second > b.second;
    }
  };
  
  auto t2 = std::chrono::high_resolution_clock::now();
  
  //DIjkstra
  
  
  //std::vector<double> Distances(NbNodes, std::numeric_limits<double>::max());                   
  
  auto t3 = std::chrono::high_resolution_clock::now();
  
  Distances[dep] = 0.0;                                                     
  
  auto t4 = std::chrono::high_resolution_clock::now();
  
  priority_queue<std::pair<int, double>, vector<std::pair<int, double> >, comp > Q;
  
  auto t5 = std::chrono::high_resolution_clock::now();
  
  Q.push(std::make_pair(dep, 0.0));   
  
  auto t6 = std::chrono::high_resolution_clock::now();
  
  //int count=0;
  std::vector<int> visited;
  //std::vector<int>::iterator it;
  while (!Q.empty()) {                                                        
    int v = Q.top().first;                                                     
    double w = Q.top().second;                                                     
    Q.pop();
    visited.push_back(v);
    if (v==node){ //node to be contracted, 
      continue;
    }
    if (Distances[v] > maxlim){
      break;
    }
    
    
    if (w <= Distances[v]) {                                                    
      
      for (unsigned int i=0; i< Graph[v].size(); i++) {
        std::pair<int,double> j = Graph[v][i];                                               
        int v2 = j.first;                                                      
        double w2 = j.second;
        if (v2==node){ //node to be contracted, 
          continue;
        }
        
        if (Distances[v] + w2 < Distances[v2]) {                               
          Distances[v2] = Distances[v] + w2;                                   
          
          Q.push(make_pair(v2, Distances[v2]));
          visited.push_back(v2);
        }
      }
    }
    
  }
  
  int total=0;
  for (int i=0; i < arr.size(); i++){
    if (Distances[arr[i]]> lim[i]) total+=1;
  }
  
  //Reinitialize distances vector
  for (int i=0; i < visited.size();i++) Distances[visited[i]]=std::numeric_limits<double>::max();
  
  return total;
  
}







int Edge_dif(int node,std::vector<std::vector<std::pair<int, double> > > &Graph,std::vector<std::vector<std::pair<int, double> > > &Graphr,int NbNodes,std::vector<double> &Distances){
  //Node degree
  int outcoming=Graph[node].size();
  int incoming=Graphr[node].size();
  
  if (incoming==0 || outcoming==0){
    return 0-incoming-outcoming;
  }
  else {
    
    
    
    //incoming > outcoming
    if (incoming <= outcoming){
      int shortcuts=0;
      for (int i=0; i < incoming; i++){
        int dep=Graphr[node][i].first;
        std::vector<double> lim(outcoming,0.0);
        double distance_1=Graphr[node][i].second;
        
        std::vector<int> arr(outcoming,0);
        for (int j=0; j < outcoming; j++){
          double distance_2=Graph[node][j].second;
          arr[j]=Graph[node][j].first;
          lim[j]=distance_1+distance_2;
          
        }
     
        shortcuts+=Dijkstra_bool(Graph,dep,arr,lim,NbNodes,node,Distances);

      }
      
      return shortcuts-incoming-outcoming;
      
    }
    //outcoming < incoming
    else {
      int shortcuts=0;
      for (int i=0; i < outcoming; i++){
        int dep=Graph[node][i].first;
        std::vector<double> lim(incoming,0.0);
        double distance_1=Graph[node][i].second;
        std::vector<int> arr(incoming,0);
        for (int j=0; j < incoming; j++){
          double distance_2=Graphr[node][j].second;
          arr[j]=Graphr[node][j].first;
          lim[j]=distance_1+distance_2;
          
        }
        
        shortcuts+=Dijkstra_bool(Graphr,dep,arr,lim,NbNodes,node,Distances);
        
      }
      
      return shortcuts-incoming-outcoming;
    }
    
    
  }
  
}


//




std::vector<std::pair<int,std::pair<int,double> > >  Dijkstra(std::vector<std::vector<std::pair<int, double> > > &Graph,//modifiable
                                                              std::vector<std::vector<std::pair<int, double> > > &Graphr,//modifiable
                                                              int dep,
                                                              std::vector<int> &arr,
                                                              std::vector<double> &lim,
                                                              int NbNodes,
                                                              int node,
                                                              std::vector<double> &Distances,
                                                              std::vector<int> &Contracted,
                                                              bool reversed,
                                                              bool &err){
  
  
  double maxlim=*max_element(lim.begin(), lim.end());
  
  
  
  auto t1 = std::chrono::high_resolution_clock::now();
  //Comparator 
  struct comp{
    
    bool operator()(const std::pair<int, double> &a, const std::pair<int, double> &b){
      return a.second > b.second;
    }
  };
  
  
  Distances[dep] = 0.0;                                                     
  
  
  priority_queue<std::pair<int, double>, vector<std::pair<int, double> >, comp > Q;
  
  Q.push(std::make_pair(dep, 0.0));   
  
  //int count=0;
  std::vector<int> visited;
  //std::vector<int>::iterator it;
  while (!Q.empty()) {                                                        
    int v = Q.top().first;                                                     
    double w = Q.top().second;                                                     
    Q.pop();
    visited.push_back(v);
    if (v==node){ //node to be contracted, 
      continue;
    }
    if (Distances[v] > maxlim){
      break;
    }
    
    
    if (w <= Distances[v]) {                                                    
      
      for (unsigned int i=0; i< Graph[v].size(); i++) {
        std::pair<int,double> j = Graph[v][i];                                               
        int v2 = j.first;                                                      
        double w2 = j.second;
        if (v2==node){ //node to be contracted, 
          continue;
        }
        
        if (Distances[v] + w2 < Distances[v2]) {                               
          Distances[v2] = Distances[v] + w2;                                   
          
          Q.push(make_pair(v2, Distances[v2]));
          visited.push_back(v2);
        }
      }
    }
    
  }
  
  
  std::vector<double> testdist;
  std::vector<std::pair<int,std::pair<int,double> > > Short;
  for (int i=0; i < arr.size(); i++){
    if (Distances[arr[i]] > lim[i]){
      if (dep==arr[i]) continue;
      
      if (reversed==true){
        
        testdist.push_back(Distances[arr[i]]);
        Short.push_back(std::make_pair(arr[i], std::make_pair(dep,lim[i])));
   
      }
      else {
        
        testdist.push_back(Distances[arr[i]]);
        Short.push_back(std::make_pair(dep, std::make_pair(arr[i],lim[i])));
        
      }
      
      
      
    }
    
  }
  
  
  
  
  //Reinitialize distances vector
  for (int i=0; i < visited.size();i++) Distances[visited[i]]=std::numeric_limits<double>::max();
  
  return Short;
  
}






void Contract(int node,
              std::vector<std::vector<std::pair<int, double> > > &Graph,//modifiable
              std::vector<std::vector<std::pair<int, double> > > &Graphr,//modifiable
              std::vector<std::vector<std::pair<int, double> > > &OrGraph,//augmented (original + shortcuts)
              
              int NbNodes,
              std::vector<double> &Distances,
              std::vector<int> &Contracted,
              int count,
              bool &err,
              std::vector<int> &ShortF,
              std::vector<int> &ShortT,
              std::vector<int> &ShortC){
  //Remove contracted nodes
  
  auto t1 = std::chrono::high_resolution_clock::now();
  
  
  
  //Node degree
  
  int outcoming=Graph[node].size();
  int incoming=Graphr[node].size();
  
  std::vector<std::pair<int,std::pair<int,double> > > Final;
  if (incoming==0 || outcoming==0){
    Contracted[node]=1;
    
  }
  else {
    
    
    
    //incoming > outcoming
    if (incoming <= outcoming){
      
      
      int shortcuts=0;
      for (int i=0; i < incoming; i++){
        
        
        int dep=Graphr[node][i].first;
        std::vector<double> lim(outcoming,0.0);
        double distance_1=Graphr[node][i].second;
        
        std::vector<int> arr(outcoming,0);
        for (int j=0; j < outcoming; j++){
          double distance_2=Graph[node][j].second;
          arr[j]=Graph[node][j].first;
          lim[j]=distance_1+distance_2;
          
        }
        
        std::vector<std::pair<int,std::pair<int,double> > > result=Dijkstra(Graph,Graphr,dep,arr,lim,NbNodes,node,Distances,Contracted,false,err);
        for (int j=0; j < result.size();j++) Final.push_back(result[j]);
        
      }
      
      
      
    }
    //outcoming < incoming
    else {
      
      
      int shortcuts=0;
      for (int i=0; i < outcoming; i++){
        int dep=Graph[node][i].first;
        std::vector<double> lim(incoming,0.0);
        double distance_1=Graph[node][i].second;
        std::vector<int> arr(incoming,0);
        for (int j=0; j < incoming; j++){
          double distance_2=Graphr[node][j].second;
          arr[j]=Graphr[node][j].first;
          lim[j]=distance_1+distance_2;
          
        }
        
        
        std::vector<std::pair<int,std::pair<int,double> > > result=Dijkstra(Graphr,Graph,dep,arr,lim,NbNodes,node,Distances,Contracted,true,err);
        for (int j=0; j < result.size();j++) Final.push_back(result[j]);
      }
      
      
    }
    
    
    
  }
  
  //modify graphs
  //Remove contracted nodes
  
  Contracted[node]=1;
  if (true) {
    
    
    for (int j=0; j < incoming; j++){
      std::vector<int> ind;
      int nd=Graphr[node][j].first;
      
      for (int i=0; i<Graph[nd].size();i++){
        if (Graph[nd][i].first==node) ind.push_back(i);
      }
      if (ind.size()>0) std::reverse(ind.begin(),ind.end());
      for (int i=0; i < ind.size(); i++) quickDelete3(ind[i], Graph[nd]);
      
      
    }
    for (int j=0; j < outcoming; j++){
      std::vector<int> ind2;
      int nd=Graph[node][j].first;
      
      for (int i=0; i<Graphr[nd].size();i++){
        if (Graphr[nd][i].first==node) ind2.push_back(i);
      }
      if (ind2.size()>0) std::reverse(ind2.begin(),ind2.end());
      for (int i=0; i < ind2.size(); i++) quickDelete3(ind2[i], Graphr[nd]);
      
      
    }
    
    
    Graph[node].erase (Graph[node].begin(),Graph[node].end());
    Graphr[node].erase (Graphr[node].begin(),Graphr[node].end());
    
  }
  
  
  
  //Shortcuts
  
  
  
  
  
  bool dup;
  for (int i=0; i<Final.size();i++){
    //Add shortcut to the final graph
    OrGraph[Final[i].first].push_back(std::make_pair(Final[i].second.first,Final[i].second.second));
    
    ShortF.push_back(Final[i].first);
    ShortT.push_back(Final[i].second.first);
    ShortC.push_back(node);
    
    
    //Duplicated
    
    dup=false;
    for (int k=0; k<Graph[Final[i].first].size(); k++){
      if (Graph[Final[i].first][k].first==Final[i].second.first && Graph[Final[i].first][k].second > Final[i].second.second){
        Graph[Final[i].first][k].second=Final[i].second.second;
        dup=true;
        break;
        
      }
    }
    
    if (dup==false) Graph[Final[i].first].push_back(std::make_pair(Final[i].second.first,Final[i].second.second));
    
    dup=false;
    for (int k=0; k<Graphr[Final[i].second.first].size(); k++){
      if (Graphr[Final[i].second.first][k].first==Final[i].first && Graphr[Final[i].second.first][k].second >Final[i].second.second){
        Graphr[Final[i].second.first][k].second=Final[i].second.second;
        dup=true;
        break;
        
      }
    }
    
    
    if (dup==false) Graphr[Final[i].second.first].push_back(std::make_pair(Final[i].first,Final[i].second.second));

    
    
  }
  
  
  
  
}
