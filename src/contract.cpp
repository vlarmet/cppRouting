#include "cgraph.h"
#include "graph.h"
#include <progress.hpp>
#include <progress_bar.hpp>
#include <vector>
#include <Rcpp.h>
using namespace std;


void quickDelete3( int idx ,std::vector<std::pair<int,double> > &vec)
{
  vec[idx] = vec.back();
  vec.pop_back();
}


void CGraph::clean(G &Graph,
                   G &OrGraph,
                   int dep, IVec &arr, DVec &lim,
                   DVec &distances){


  double maxlim=*max_element(lim.begin(), lim.end());


  //Dijkstra

  distances[dep] = 0.0;
  PQ Q;
  Q.push(make_pair(dep, 0.0));
  IVec visited;

  set<int> visited2;
  set<int> sel;
  for (int i = 0; i < arr.size(); i++) sel.insert(arr[i]);

  while (!Q.empty()) {
    int v = Q.top().first;
    double w = Q.top().second;
    Q.pop();
    visited.push_back(v);

    if (sel.find(v) != sel.end()) visited2.insert(v);

    if (visited2.size() == arr.size()) break;

    if (distances[v] > maxlim){
      break;
    }


    if (w <= distances[v]) {

      for (unsigned int i=0; i< OrGraph[v].size(); i++) {
        pair<int,double> j = OrGraph[v][i];
        int v2 = j.first;
        double w2 = j.second;

        if (distances[v] + w2 < distances[v2]) {
          distances[v2] = distances[v] + w2;


          Q.push(make_pair(v2, distances[v2]));
          visited.push_back(v2);
        }
      }
    }

  }


  for (int i=0; i < arr.size(); i++){

    if (distances[arr[i]]< lim[i]) {

      for (int j=0; j < Graph[dep].size();j++){
        if (Graph[dep][j].first==arr[i]){
          quickDelete3(j, Graph[dep]);
          break;
        }
      }

    }
  }

  //Reinitialize distances vector
  for (int i=0; i < visited.size();i++) distances[visited[i]]=numeric_limits<double>::max();

}



int CGraph::find_shortcuts(G &graph,
                           int dep, IVec &arr,
                           DVec &lim, int node,
                           DVec &distances){
  double maxlim=*max_element(lim.begin(), lim.end());

  distances[dep] = 0.0;
  PQ Q;
  Q.push(make_pair(dep, 0.0));
  IVec visited;

  while (!Q.empty()) {
    int v = Q.top().first;
    double w = Q.top().second;
    Q.pop();
    visited.push_back(v);
    if (v==node){ //node to be contracted,
      continue;
    }
    if (distances[v] > maxlim){
      break;
    }

    if (w <= distances[v]) {

      for (unsigned int i=0; i< graph[v].size(); i++) {
        pair<int,double> j = graph[v][i];
        int v2 = j.first;
        double w2 = j.second;
        if (v2==node){ //node to be contracted,
          continue;
        }

        if (distances[v] + w2 < distances[v2]) {
          distances[v2] = distances[v] + w2;

          Q.push(make_pair(v2, distances[v2]));
          visited.push_back(v2);
        }
      }
    }

  }

  int total = 0;
  for (int i=0; i < arr.size(); i++){
    if (distances[arr[i]]> lim[i]) total+=1;
  }

  //Reinitialize distances vector
  for (int i=0; i < visited.size();i++) distances[visited[i]] = numeric_limits<double>::max();

  return total;
}


int CGraph::edge_dif(int node, G &graph, G &graphr,
                     DVec &distances){
  //Node degree

  int outcoming = graph[node].size();
  int incoming = graphr[node].size();

  if (incoming==0 || outcoming==0){
    return 0-incoming-outcoming;
  }
  else {



    //incoming > outcoming
    if (incoming <= outcoming){
      int shortcuts=0;
      for (int i=0; i < incoming; i++){
        int dep = graphr[node][i].first;
        DVec lim(outcoming,0.0);
        double distance_1 = graphr[node][i].second;

        IVec arr(outcoming,0);
        for (int j=0; j < outcoming; j++){
          double distance_2 = graph[node][j].second;
          arr[j] = graph[node][j].first;
          lim[j] = distance_1+distance_2;

        }

        shortcuts += find_shortcuts(graph, dep, arr, lim, node, distances);

      }

      return shortcuts-incoming-outcoming;

    }

    //outcoming < incoming
    else {
      int shortcuts=0;
      for (int i=0; i < outcoming; i++){
        int dep = graph[node][i].first;
        DVec lim(incoming,0.0);
        double distance_1 = graph[node][i].second;
        IVec arr(incoming,0);
        for (int j=0; j < incoming; j++){
          double distance_2 = graphr[node][j].second;
          arr[j] = graphr[node][j].first;
          lim[j] = distance_1+distance_2;

        }

        shortcuts += find_shortcuts(graphr, dep, arr, lim, node, distances);

      }

      return shortcuts-incoming-outcoming;
    }


  }

}


vector<pair<int, pair<int,double>>> CGraph::get_shortcuts(G &graph,
                                                          int dep,
                                                          IVec &arr,
                                                          DVec &lim,
                                                          int node,
                                                          DVec &distances,
                                                          bool reversed){
  double maxlim=*max_element(lim.begin(), lim.end());

  distances[dep] = 0.0;
  PQ Q;
  Q.push(make_pair(dep, 0.0));
  IVec visited;

  while (!Q.empty()) {
    int v = Q.top().first;
    double w = Q.top().second;
    Q.pop();
    visited.push_back(v);
    if (v==node){ //node to be contracted,
      continue;
    }
    if (distances[v] > maxlim){
      break;
    }


    if (w <= distances[v]) {

      for (unsigned int i=0; i< graph[v].size(); i++) {
        pair<int,double> j = graph[v][i];
        int v2 = j.first;
        double w2 = j.second;
        if (v2==node){ //node to be contracted,
          continue;
        }

        if (distances[v] + w2 < distances[v2]) {
          distances[v2] = distances[v] + w2;

          Q.push(make_pair(v2, distances[v2]));
          visited.push_back(v2);
        }
      }
    }

  }

  DVec testdist;
  vector<pair<int,pair<int,double> > > Short;
  for (int i=0; i < arr.size(); i++){
    if (distances[arr[i]] > lim[i]){
      if (dep==arr[i]) continue;

      if (reversed==true){

        testdist.push_back(distances[arr[i]]);
        Short.push_back(make_pair(arr[i], make_pair(dep,lim[i])));

      }
      else {

        testdist.push_back(distances[arr[i]]);
        Short.push_back(make_pair(dep, make_pair(arr[i],lim[i])));

      }
    }
  }

  //Reinitialize distances vector
  for (int i = 0; i < visited.size(); i++) distances[visited[i]]=numeric_limits<double>::max();

  return Short;
}


void CGraph::contract_one_node(int node,
                               G &graph,
                               G &graphr,
                               G &cgraph,
                               DVec &distances,
                               IVec &contracted){

  //Node degree

  int outcoming = graph[node].size();
  int incoming = graphr[node].size();

  vector<pair<int,pair<int,double> > > Final;
  if (incoming==0 || outcoming==0){
    contracted[node] = 1;

  }
  else {



    //incoming > outcoming
    if (incoming <= outcoming){


      int shortcuts=0;
      for (int i=0; i < incoming; i++){


        int dep = graphr[node][i].first;
        DVec lim(outcoming,0.0);
        double distance_1 = graphr[node][i].second;

        IVec arr(outcoming,0);
        for (int j=0; j < outcoming; j++){
          double distance_2 = graph[node][j].second;
          arr[j] = graph[node][j].first;
          lim[j] = distance_1+distance_2;

        }

        vector<pair<int,pair<int,double> > > result = get_shortcuts(graph, dep, arr, lim, node, distances,false);
        for (int j=0; j < result.size();j++) Final.push_back(result[j]);

      }



    }
    //outcoming < incoming
    else {


      int shortcuts=0;
      for (int i=0; i < outcoming; i++){
        int dep = graph[node][i].first;
        DVec lim(incoming,0.0);
        double distance_1 = graph[node][i].second;
        IVec arr(incoming,0);
        for (int j=0; j < incoming; j++){
          double distance_2 = graphr[node][j].second;
          arr[j] = graphr[node][j].first;
          lim[j] = distance_1+distance_2;

        }


        vector<pair<int,pair<int,double> > > result = get_shortcuts(graphr, dep, arr, lim, node, distances,true);
        for (int j=0; j < result.size();j++) Final.push_back(result[j]);
      }

    }

  }

  //modify graphs
  //Remove contracted nodes

  contracted[node]=1;

  for (int j=0; j < incoming; j++){
    IVec ind;
    int nd = graphr[node][j].first;

    for (int i=0; i < graph[nd].size();i++){
      if (graph[nd][i].first==node) ind.push_back(i);
    }
    if (ind.size() > 0) reverse(ind.begin(),ind.end());
    for (int i=0; i < ind.size(); i++) quickDelete3(ind[i], graph[nd]);


  }
  for (int j=0; j < outcoming; j++){
    IVec ind2;
    int nd = graph[node][j].first;

    for (int i=0; i < graphr[nd].size();i++){
      if (graphr[nd][i].first==node) ind2.push_back(i);
    }
    if (ind2.size()>0) reverse(ind2.begin(),ind2.end());
    for (int i=0; i < ind2.size(); i++) quickDelete3(ind2[i], graphr[nd]);


  }


  graph[node].erase (graph[node].begin(),graph[node].end());
  graphr[node].erase (graphr[node].begin(),graphr[node].end());




  //Shortcuts

  bool dup;
  for (int i=0; i < Final.size(); i++){
    //Add shortcut to the final graph
    cgraph[Final[i].first].push_back(make_pair(Final[i].second.first,Final[i].second.second));

    shortf.push_back(Final[i].first);
    shortt.push_back(Final[i].second.first);
    shortc.push_back(node);


    //Duplicated

    dup=false;

    for (int k=0; k < graph[Final[i].first].size(); k++){
      if (graph[Final[i].first][k].first==Final[i].second.first && graph[Final[i].first][k].second > Final[i].second.second){
        graph[Final[i].first][k].second=Final[i].second.second;
        dup=true;
        break;

      }
    }

    if (dup==false) graph[Final[i].first].push_back(make_pair(Final[i].second.first,Final[i].second.second));

    dup=false;
    for (int k=0; k < graphr[Final[i].second.first].size(); k++){
      if (graphr[Final[i].second.first][k].first==Final[i].first && graphr[Final[i].second.first][k].second >Final[i].second.second){
        graphr[Final[i].second.first][k].second=Final[i].second.second;
        dup=true;
        break;

      }
    }


    if (dup==false) graphr[Final[i].second.first].push_back(make_pair(Final[i].first,Final[i].second.second));

  }
}

// [[Rcpp::depends(RcppProgress)]]
void CGraph::contract(bool display_progress){

  DVec distances(nbnode, numeric_limits<double>::max());
  G data = cdata; // modifiable graph

  if (display_progress) Rcpp::Rcout << "Cleaning graph..." <<endl;


  for (int i = 0; i< nbnode; i++){
    IVec arr;
    DVec lim;
    for (int j=0; j < cdata[i].size(); j++){
      arr.push_back(cdata[i][j].first);
      lim.push_back(cdata[i][j].second);

    }
    if (arr.size()==0) continue;

    clean(cdata,
          data,
          i,
          arr,
          lim,
          distances);
  }


  //Reversed graph
  G dataR = getReverse(cdata);
  data = cdata;


  //ordering
  if (display_progress) Rcpp::Rcout << "Initial node ordering..."<< endl;
  priority_queue<pair<int, int>, vector<pair<int, int> >, comp > Queue;

  IVec contracted(nbnode, 0);



  for (int i=0; i < nbnode; i++){
    Rcpp::checkUserInterrupt();


    int ED = edge_dif(i, cdata, dataR, distances);



    Queue.push(make_pair(i, ED));
  }

  //Contracting


  int count=0;
  int count2=0;
  int count3=0;
  bool error=false;

  if (display_progress) Rcpp::Rcout << "Contracting nodes..."<< endl;
  Progress p(nbnode, display_progress);


  //Contracting all nodes
  while(!Queue.empty()){

    if (error) break;

    Rcpp::checkUserInterrupt();

    int v=Queue.top().first;
    int imp=Queue.top().second;
    Queue.pop();

    if (contracted[v] == 1) {
      count3+=1;
      continue; //already contracted
    }


    int DN = 0;
    for (int i = 0; i < cdata[v].size(); i++){
      if (contracted[cdata[v][i].first] ==1 ) DN += 1;
    }


    int New_imp = edge_dif(v,data, dataR, distances) * 190 + DN * 120;


    if (Queue.empty()){
      count+=1;
      //   p.increment();

      rank[v]=count;
      contracted[v]=1;

      contract_one_node(v,
                        data,
                        dataR,
                        cdata,
                        distances,
                        contracted);
    }

    else{
      if (New_imp > Queue.top().second){ //compare to the next min
        count2+=1;
        //reinsert
        Queue.push(make_pair(v,New_imp));

      }


      else{

        //contract

        count+=1;
        p.increment();
        rank[v] = count;
        contracted[v] = 1;


        contract_one_node(v,
                          data,
                          dataR,
                          cdata,
                          distances,
                          contracted);


      }
    }
  }

  // update number of edges
  count=0;
  for (int i = 0; i < nbnode; i++) count += cdata[i].size();
  nbedge = count;

}

