#include "graph.h"
#include "cgraph.h"
#include "aon.h"
#include "stall.h"
#include <string>
#include <RcppParallel.h>
#include <Rcpp.h>
using namespace std;
using namespace RcppParallel;

aonGraph::aonGraph(Graph* gr, IVec dep, IVec arr, DVec demand, int algo) :
  m_gr(gr), m_dep(dep), m_arr(arr), m_demand(demand), algorithm(algo)
{
  m_result.resize(m_gr->nbedge, 0.0);

  if (algorithm == 0 || algorithm == 1){
    // first, a map
    map<int, pair<IVec, DVec>> od_map;
    map<int, pair<IVec, DVec>>::iterator ite;

    for (int i = 0; i < m_dep.size(); i++) {
      if (algorithm == 0){
        od_map[m_dep[i]].first.push_back(m_arr[i]);
        od_map[m_dep[i]].second.push_back(m_demand[i]);
      }
      if (algorithm == 1){
        od_map[m_arr[i]].first.push_back(m_dep[i]);
        od_map[m_arr[i]].second.push_back(m_demand[i]);
      }

    }

    // then, convert it to a vector

    vector<pair<int, pair<vector<int>, vector<double>>>> od(od_map.size());
    int ind = 0;
    for (ite = od_map.begin(); ite != od_map.end(); ite++){
      od[ind] = make_pair(ite->first, ite->second);
      ind++;
    }

    m_od = od;
    n_iter = od.size();
  }
  if (algorithm == 2 || algorithm == 3){
    m_od.resize(0);
    n_iter = m_dep.size();
  }
}

aonGraph::aonGraph(aonGraph &pn, RcppParallel::Split):
  m_gr(pn.m_gr), m_dep(pn.m_dep), m_arr(pn.m_arr), m_demand(pn.m_demand), algorithm(pn.algorithm),
  m_result(pn.m_result), m_od(pn.m_od), n_iter(pn.n_iter)
{

}

void aonGraph::operator()(std::size_t begin, std::size_t end){
  if (algorithm == 0) dijkstra(begin, end);
  if (algorithm == 1) dijkstra_reverse(begin, end);
  if (algorithm == 2) bidir(begin, end);
  if (algorithm == 3) nba(begin, end);
}

void aonGraph::join(aonGraph &pb){
  for (int i = 0; i < m_result.size(); i++){
    m_result[i] += pb.m_result[i];
  }
}

void aonGraph::dijkstra(std::size_t begin, std::size_t end){
  std::vector<double> Distances(m_gr->nbnode, std::numeric_limits<double>::max());
  std::vector<int> Parents(m_gr->nbnode, -1);


  for (std::size_t k=begin; k!=end;k++){



    int StartNode=m_od[k].first;

    Distances[StartNode] = 0.0;
    PQ Q;
    Q.push(std::make_pair(StartNode, 0.0));


    while (!Q.empty()) {
      int v = Q.top().first;
      double w = Q.top().second;
      Q.pop();

      if (w <= Distances[v]) {

        for (int i=m_gr->indG[v]; i< m_gr->indG[v+1]; i++){
          int v2 = m_gr->nodeG[i];
          double w2 = m_gr->wG[i];

          if (Distances[v] + w2 < Distances[v2]) {
            Distances[v2] = Distances[v] + w2;
            Parents[v2] = v;
            Q.push(std::make_pair(v2, Distances[v2]));
          }
        }
      }

    }

    IVec arr = m_od[k].second.first;
    DVec dem = m_od[k].second.second;

    // load flow
    for (int i = 0; i < arr.size(); i++){
      int node = arr[i];
      if (Distances[node] == std::numeric_limits<double>::max()) continue;

      for (auto p = node; p != -1; p = Parents[p]){
        int prev = Parents[p];
        if (prev == -1) break;

        int edge_index = -1;
        double w = numeric_limits<double>::max();

        for (int j = m_gr->indG[prev]; j < m_gr->indG[prev + 1]; j++){
          if (m_gr->nodeG[j] == p){
            if (m_gr->wG[j] < w){
              edge_index = j;
              w = m_gr->wG[j];
            }

          }
        }
        if (edge_index < 0) Rcpp::Rcout << "ok" << endl;
        m_result[edge_index] += dem[i];
      }
    }

    std::fill(Distances.begin(), Distances.end(), std::numeric_limits<double>::max());
    std::fill(Parents.begin(), Parents.end(), -1);

  }
}

void aonGraph::dijkstra_reverse(std::size_t begin, std::size_t end){
  std::vector<double> Distances(m_gr->nbnode, std::numeric_limits<double>::max());
  std::vector<int> Parents(m_gr->nbnode, -1);


  for (std::size_t k=begin; k!=end;k++){



    int StartNode=m_od[k].first;

    Distances[StartNode] = 0.0;
    PQ Q;
    Q.push(std::make_pair(StartNode, 0.0));


    while (!Q.empty()) {
      int v = Q.top().first;
      double w = Q.top().second;
      Q.pop();

      if (w <= Distances[v]) {

        for (int i=m_gr->indGr[v]; i< m_gr->indGr[v+1]; i++){
          int v2 = m_gr->nodeGr[i];
          double w2 = m_gr->wGr[i];

          if (Distances[v] + w2 < Distances[v2]) {
            Distances[v2] = Distances[v] + w2;
            Parents[v2] = v;
            Q.push(std::make_pair(v2, Distances[v2]));
          }
        }
      }

    }

    IVec arr = m_od[k].second.first;
    DVec dem = m_od[k].second.second;

    // load flow
    for (int i = 0; i < arr.size(); i++){
      int node = arr[i];
      if (Distances[node] == std::numeric_limits<double>::max()) continue;

      for (auto p = node; p != -1; p = Parents[p]){
        int prev = Parents[p];
        if (prev == -1) break;

        int edge_index = -1;
        double w = numeric_limits<double>::max();

        for (int j = m_gr->indG[p]; j < m_gr->indG[p + 1]; j++){
          if (m_gr->nodeG[j] == prev){
            if (m_gr->wG[j] < w){
              edge_index = j;
              w = m_gr->wG[j];
            }

          }
        }
        m_result[edge_index] += dem[i];
      }
    }

    std::fill(Distances.begin(), Distances.end(), std::numeric_limits<double>::max());
    std::fill(Parents.begin(), Parents.end(), -1);

  }
}

void aonGraph::bidir(std::size_t begin, std::size_t end){
  std::vector<double> Distances(m_gr->nbnode, std::numeric_limits<double>::max());
  std::vector<double> Distances2(m_gr->nbnode, std::numeric_limits<double>::max());
  std::vector<int> Parents(m_gr->nbnode, -1);
  std::vector<int> Parents2(m_gr->nbnode, -1);
  std::vector <int> Visited(m_gr->nbnode, 0);
  std::vector <int> Visiting(m_gr->nbnode, 0);
  std::vector <int> Visited2(m_gr->nbnode, 0);
  std::vector <int> Visiting2(m_gr->nbnode, 0);

  for (std::size_t k=begin; k!=end;k++){



    int StartNode=m_dep[k];
    int EndNode=m_arr[k];
    Distances[StartNode] = 0.0;
    Distances2[EndNode] = 0.0;

    PQ Q;
    PQ Qr;
    Q.push(std::make_pair(StartNode, 0.0));
    Qr.push(std::make_pair(EndNode, 0.0));
    Visiting[StartNode]=1;
    Visiting2[EndNode]=1;

    int mid = -1;
    double mu=std::numeric_limits<double>::max();


    while (!Q.empty() && !Qr.empty()) {
      if (Q.top().second+Qr.top().second >= mu){
        break;
      }


      if (!Q.empty()){
        int v=Q.top().first;
        int w=Q.top().second;
        Q.pop();


        Visited[v]=1;


        if (w <= Distances[v]) {
          for (int i=m_gr->indG[v]; i< m_gr->indG[v+1]; i++){
            int v2 = m_gr->nodeG[i];
            double w2 = m_gr->wG[i];


            if (Distances[v] + w2 < Distances[v2]) {
              Distances[v2] = Distances[v] + w2;
              Parents[v2] = v;
              Q.push(std::make_pair(v2, Distances[v2]));
              Visiting[v2]=1;

            }
          }
        }
        if ((Visited2[v]==1 || Visiting2[v]==1)  && (Distances[v]+Distances2[v]) < mu){
          mid = v;
          mu=Distances[v]+Distances2[v];

        }
      }

      if (!Qr.empty()){
        int vv=Qr.top().first;
        int ww=Qr.top().second;
        Qr.pop();


        Visited2[vv]=1;

        if (ww <= Distances2[vv]) {
          for (int i=m_gr->indGr[vv]; i< m_gr->indGr[vv+1]; i++){
            int vv2 = m_gr->nodeGr[i];
            double ww2 = m_gr->wGr[i];


            if (Distances2[vv] + ww2 < Distances2[vv2]) {
              Distances2[vv2] = Distances2[vv] + ww2;
              Parents2[vv2] = vv;
              Qr.push(std::make_pair(vv2, Distances2[vv2]));
              Visiting2[vv2]=1;
            }
          }
        }
        if ((Visited[vv]==1 || Visiting[vv]==1) && (Distances[vv]+Distances2[vv]) < mu){
          mid = vv;
          mu=Distances[vv]+Distances2[vv];

        }
      }

    }

    if (mid != -1) {
      // path

      IVec result2;
      for (auto p = Parents2[mid]; p != -1; p = Parents2[p]){
        result2.insert(result2.begin(),p);
      }

      if (Distances[mid]!=numeric_limits<double>::max() || Distances2[mid]!=numeric_limits<double>::max()){
        result2.push_back(mid);
      }

      for (auto p = Parents[mid]; p != -1; p = Parents[p]){
        result2.push_back(p);
      }

      std::reverse(result2.begin(), result2.end());

      // update aon
      if (result2.size() > 0){

        int pos = 0;
        int dep = result2[pos];
        int second = result2[pos + 1];
        int arr = result2.back();


        while (dep != arr){

          second = result2[pos + 1];
          int edge_index = -1;
          double w = std::numeric_limits<double>::max();
          for (int i = m_gr->indG[dep]; i < m_gr->indG[dep + 1]; i++){

            if (m_gr->nodeG[i] == second && m_gr->wG[i] < w){
              w =  m_gr->wG[i];
              edge_index = i;

            }
          }
          m_result[edge_index] += m_demand[k];
          dep = second;
          pos += 1;
          //break;
        }
      }
    }

    std::fill(Distances.begin(), Distances.end(), std::numeric_limits<double>::max());
    std::fill(Distances2.begin(), Distances2.end(), std::numeric_limits<double>::max());
    std::fill(Parents.begin(), Parents.end(), -1);
    std::fill(Parents2.begin(), Parents2.end(), -1);
    std::fill(Visited.begin(), Visited.end(), 0);
    std::fill(Visited2.begin(), Visited2.end(), 0);
    std::fill(Visiting.begin(), Visiting.end(), 0);
    std::fill(Visiting2.begin(), Visiting2.end(), 0);

  }
}

void aonGraph::nba(std::size_t begin, std::size_t end){
  std::vector<double> Distances(m_gr->nbnode, std::numeric_limits<double>::max());
  std::vector<double> Distances2(m_gr->nbnode, std::numeric_limits<double>::max());
  std::vector<int> Parents(m_gr->nbnode, -1);
  std::vector<int> Parents2(m_gr->nbnode, -1);
  std::vector <int> Visited(m_gr->nbnode,0);
  std::vector <int> Visited1check(m_gr->nbnode,0);
  std::vector <int> Visited2check(m_gr->nbnode,0);

  for (std::size_t k=begin; k!=end;k++){



    int StartNode=m_dep[k];
    int EndNode=m_arr[k];

    if (StartNode==EndNode){
      continue;
    }

    double lata=m_gr->lat[EndNode];
    double lona=m_gr->lon[EndNode];
    double lata2=m_gr->lat[StartNode];
    double lona2=m_gr->lon[StartNode];

    Distances[StartNode] = 0.0;
    Visited1check[StartNode]=1;
    Distances2[EndNode] = 0.0;
    Visited2check[EndNode]=1;
    PQ Q;
    PQ Qr;
    Q.push(std::make_pair(StartNode, sqrt(pow(m_gr->lat[StartNode]-lata,2)+pow(m_gr->lon[StartNode]-lona,2))/m_gr->k));
    Qr.push(std::make_pair(EndNode, sqrt(pow(m_gr->lat[EndNode]-lata2,2)+pow(m_gr->lon[EndNode]-lona2,2))/m_gr->k));

    double total1=sqrt(pow(lata2-lata,2)+pow(lona2-lona,2))/m_gr->k;
    double total2=sqrt(pow(lata2-lata,2)+pow(lona2-lona,2))/m_gr->k;
    int mid = -1;
    double mu=std::numeric_limits<double>::max();

    while (!Q.empty() && !Qr.empty()) {
      //Forward
      if (Q.size() < Qr.size()){
        int v=Q.top().first;
        Q.pop();
        if (Visited[v]==0){
          Visited[v]=1;

          if ((Distances[v] + sqrt(pow(m_gr->lat[v]-lata,2)+pow(m_gr->lon[v]-lona,2))/m_gr->k) >= mu || (Distances[v] + total2 - sqrt(pow(m_gr->lat[v]-lata2,2)+pow(m_gr->lon[v]-lona2,2))/m_gr->k) >= mu){}

          else {
            for (int i=m_gr->indG[v]; i< m_gr->indG[v+1]; i++){
              int v2 = m_gr->nodeG[i];
              double w2 = m_gr->wG[i];
              if (Visited[v2]==1){
                continue;
              }
              double tentative=Distances[v]+w2;

              if (Visited1check[v2]==0  || Distances[v2] > tentative){
                Distances[v2]=tentative;
                Visited1check[v2]=1;
                Parents[v2]=v;
                Q.push(std::make_pair(v2, tentative + sqrt(pow(m_gr->lat[v2]-lata,2)+pow(m_gr->lon[v2]-lona,2))/m_gr->k));

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

          if ((Distances2[vv] + sqrt(pow(m_gr->lat[vv]-lata2,2)+pow(m_gr->lon[vv]-lona2,2))/m_gr->k) >= mu || (Distances2[vv] + total1 - sqrt(pow(m_gr->lat[vv]-lata,2)+pow(m_gr->lon[vv]-lona,2))/m_gr->k) >= mu){}

          else {
            for (int i=m_gr->indGr[vv]; i< m_gr->indGr[vv+1]; i++){
              int vv2 = m_gr->nodeGr[i];
              double ww2 = m_gr->wGr[i];
              if (Visited[vv2]==1){
                continue;
              }
              double tentative=Distances2[vv]+ww2;

              if (Visited2check[vv2]==0  || Distances2[vv2] > tentative){
                Distances2[vv2]=tentative;
                Visited2check[vv2]=1;
                Parents2[vv2]=vv;
                Qr.push(std::make_pair(vv2, tentative + sqrt(pow(m_gr->lat[vv2]-lata2,2)+pow(m_gr->lon[vv2]-lona2,2))/m_gr->k));

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

    if (mid != -1) {
      // path
      IVec result2;
      for (auto p = Parents2[mid]; p != -1; p = Parents2[p]){
        result2.insert(result2.begin(),p);
      }

      if (Distances[mid]!=numeric_limits<double>::max() || Distances2[mid]!=numeric_limits<double>::max()){
        result2.push_back(mid);
      }

      for (auto p = Parents[mid]; p != -1; p = Parents[p]){
        result2.push_back(p);
      }

      std::reverse(result2.begin(), result2.end());


      // update aon
      if (result2.size() > 0){

        int pos = 0;
        int dep = result2[pos];
        int second = result2[pos + 1];
        int arr = result2.back();


        while (dep != arr){

          second = result2[pos + 1];
          int edge_index = -1;
          double w = std::numeric_limits<double>::max();
          for (int i = m_gr->indG[dep]; i < m_gr->indG[dep + 1]; i++){


            if (m_gr->nodeG[i] == second && m_gr->wG[i] < w){
              w = m_gr->wG[i];
              edge_index = i;

            }
          }
          m_result[edge_index] += m_demand[k];
          dep = second;
          pos += 1;
          //break;
        }
      }
    }





    //Reinitialize
    std::fill(Distances.begin(),Distances.end(),std::numeric_limits<double>::max());
    std::fill(Distances2.begin(),Distances2.end(),std::numeric_limits<double>::max());
    std::fill(Parents.begin(),Parents.end(),-1);
    std::fill(Parents2.begin(),Parents2.end(),-1);
    std::fill(Visited.begin(),Visited.end(),0);
    std::fill(Visited1check.begin(),Visited1check.end(),0);
    std::fill(Visited2check.begin(),Visited2check.end(),0);

  }
}



/////////// contracted

aonGraphC::aonGraphC(CGraph* gr, Graph* original_graph, IVec dep, IVec arr, DVec demand, int algo) :
  m_gr(gr), m_or(original_graph), m_dep(dep), m_arr(arr), m_demand(demand), algorithm(algo)
{
  m_result.resize(m_or->nbedge, 0.0);

  if (algorithm == 0 || algorithm == 1){
    // first, a map
    map<int, pair<IVec, DVec>> od_map;
    map<int, pair<IVec, DVec>>::iterator ite;

    for (int i = 0; i < m_dep.size(); i++) {
      if (algorithm == 0){
        od_map[m_dep[i]].first.push_back(m_arr[i]);
        od_map[m_dep[i]].second.push_back(m_demand[i]);
      }
      if (algorithm == 1){
        od_map[m_arr[i]].first.push_back(m_dep[i]);
        od_map[m_arr[i]].second.push_back(m_demand[i]);
      }

    }

    // then, convert it to a vector

    vector<pair<int, pair<vector<int>, vector<double>>>> od(od_map.size());
    int ind = 0;
    for (ite = od_map.begin(); ite != od_map.end(); ite++){
      od[ind] = make_pair(ite->first, ite->second);
      ind++;
    }

    m_od = od;
    n_iter = od.size();
  }
  if (algorithm == 2){
    m_od.resize(0);
    n_iter = m_dep.size();
  }

}

aonGraphC::aonGraphC(aonGraphC &pn, RcppParallel::Split):
  m_gr(pn.m_gr), m_or(pn.m_or), m_dep(pn.m_dep), m_arr(pn.m_arr), m_demand(pn.m_demand), algorithm(pn.algorithm),
  m_result(pn.m_result), m_od(pn.m_od), n_iter(pn.n_iter)
{

}

void aonGraphC::operator()(std::size_t begin, std::size_t end){
  if (algorithm == 0) phast(begin, end);
  if (algorithm == 1) phastr(begin, end);
  if (algorithm == 2) dijkstra_bidir(begin, end);

}

void aonGraphC::join(aonGraphC &pb){
  for (int i = 0; i < m_result.size(); i++){
    m_result[i] += pb.m_result[i];
  }
}

void aonGraphC::dijkstra_bidir(std::size_t begin, std::size_t end){
  DVec distances(m_gr->nbnode, numeric_limits<double>::max());
  DVec distances2(m_gr->nbnode, numeric_limits<double>::max());
  IVec visited1(m_gr->nbnode, 0);
  IVec visited2(m_gr->nbnode, 0);
  IVec visited;
  IVec parents(m_gr->nbnode, -1);
  IVec parents2(m_gr->nbnode, -1);


  for (std::size_t k = begin; k != end; k++){



    int start = m_dep[k];
    int end = m_arr[k];
    distances[start] = 0.0;
    distances2[end] = 0.0;
    PQ Q;
    PQ Qr;
    Q.push(make_pair(start, 0.0));
    Qr.push(make_pair(end, 0.0));
    int mid;
    double mu=numeric_limits<double>::max();

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

        visited.push_back(v);

        visited1[v]=1;

        if ((visited2[v]==1)  && (distances[v]+distances2[v]) < mu){
          mid=v;
          mu=distances[v]+distances2[v];

        }

        if (w <= distances[v] && !Stall_par(v, distances, m_gr->nodeGr, m_gr->wGr, m_gr->indGr)) {
          for (int i=m_gr->indG[v]; i< m_gr->indG[v+1]; i++){
            int v2 = m_gr->nodeG[i];
            double w2 = m_gr->wG[i];

            if (distances[v] + w2 < distances[v2]) {
              distances[v2] = distances[v] + w2;
              Q.push(make_pair(v2, distances[v2]));
              visited1[v2]=1;
              parents[v2] = v;
              visited.push_back(v2);

            }
          }
        }

      }

      if (!Qr.empty()){
        int vv=Qr.top().first;
        double ww=Qr.top().second;
        Qr.pop();

        visited.push_back(vv);

        visited2[vv]=1;


        if ((visited1[vv]== 1) && (distances[vv]+distances2[vv]) < mu){
          mid=vv;
          mu=distances[vv]+distances2[vv];

        }

        if (ww <= distances2[vv] && !Stall_par(vv, distances2, m_gr->nodeG, m_gr->wG, m_gr->indG)) {
          for (int i=m_gr->indGr[vv]; i< m_gr->indGr[vv+1]; i++){
            int vv2 = m_gr->nodeGr[i];
            double ww2 = m_gr->wGr[i];

            if (distances2[vv] + ww2 < distances2[vv2]) {
              distances2[vv2] = distances2[vv] + ww2;

              Qr.push(make_pair(vv2, distances2[vv2]));
              visited2[vv2]=1;
              parents2[vv2] = vv;
              visited.push_back(vv2);

            }
          }
        }

      }

    }

    //Extracting path
    IVec result2;
    if (mu < numeric_limits<double>::max() ){

      for (auto p = parents2[mid]; p != -1; p = parents2[p]){
        result2.insert(result2.begin(),p);
      }

      if (distances[mid]!=numeric_limits<double>::max() || distances2[mid]!=numeric_limits<double>::max()){
        result2.push_back(mid);
      }

      for (auto p = parents[mid]; p != -1; p = parents[p]){
        result2.push_back(p);
      }

      reverse(result2.begin(),result2.end());


      if (result2.size() > 1){

        //m_gr->unpack(result2);

      }

      // load aon
      if (result2.size() > 1){

        int pos = 0;
        int dep = result2[pos];
        int second = result2[pos + 1];
        int arr = result2.back();


        while (dep != arr){

          second = result2[pos + 1];
          int edge_index = -1;
          double w = std::numeric_limits<double>::max();
          for (int i = m_or->indG[dep]; i < m_or->indG[dep + 1]; i++){

            if (m_or->nodeG[i] == second && m_or->wG[i] < w){
              w =  m_or->wG[i];
              edge_index = i;

            }
          }


          m_result[edge_index] += m_demand[k];
          dep = second;
          pos += 1;
        }
      }


    }


    for (int i = 0;  i <visited.size();i++){
      visited1[visited[i]]=0;
      visited2[visited[i]]=0;
      distances[visited[i]]=numeric_limits<double>::max();
      distances2[visited[i]]=numeric_limits<double>::max();
      parents[visited[i]]= -1;
      parents2[visited[i]]= -1;
    }
    visited.clear();


  }
}

void aonGraphC::phast(std::size_t begin, std::size_t end){
  DVec Distances(m_gr->nbnode, std::numeric_limits<double>::max());
  IVec parents(m_gr->nbnode, -1);
  IVec parents2(m_gr->nbnode, -1);

  for (std::size_t k=begin; k!=end;k++){



    int StartNode = m_od[k].first;

    Distances[StartNode] = 0.0;
    PQ Q;
    Q.push(std::make_pair(StartNode, 0.0));


    while (true) {

      if (Q.empty()){
        break;
      }

      if (!Q.empty()){

        int v=Q.top().first;
        double w=Q.top().second;
        Q.pop();

        if (w <= Distances[v]) {
          if (Stall_par(v, Distances, m_gr->nodeGr, m_gr->wGr, m_gr->indGr)) continue;

          for (int i = m_gr->indG[v]; i< m_gr->indG[v+1]; i++){
            int v2 = m_gr->nodeG[i];
            double w2 = m_gr->wG[i];


            if (Distances[v] + w2 < Distances[v2]) {
              Distances[v2] = Distances[v] + w2;
              parents[v2] = v;
              Q.push(std::make_pair(v2, Distances[v2]));

            }
          }
        }

      }
    }

    //Backward

    for (int i=0; i < (m_gr->indGr.size()-1); i++){

      for (int j=m_gr->indGr[i]; j < m_gr->indGr[i+1]; j++){


        if (Distances[m_gr->nodeGr[j]]+m_gr->wGr[j] < Distances[i]) {
          Distances[i] = Distances[m_gr->nodeGr[j]]+m_gr->wGr[j];
          parents2[i] = m_gr->nodeGr[j];
        }

      }
    }

    IVec arr_vec = m_od[k].second.first;
    DVec dem = m_od[k].second.second;

    for (unsigned int i=0; i < arr_vec.size();i++){
      int end = arr_vec[i];

      IVec result2;
      int last_node = -1;
      for (auto p = end; p != -1; p = parents2[p]){
        result2.push_back(p);
        last_node = p;
      }


      reverse(result2.begin(),result2.end());

      if (last_node != -1){
        for (auto p = parents[last_node]; p != -1; p = parents[p]){
          result2.insert(result2.begin(), p);
          last_node = p;
        }
      }

      if (result2.size() > 1){
    //    m_gr->unpack(result2);
      }

      // load AON
      if (result2.size() > 1){

        int pos = 0;
        int dep = result2[pos];
        int second = result2[pos + 1];
        int arr = result2.back();


        while (dep != arr){

          second = result2[pos + 1];
          int edge_index = -1;
          double w = std::numeric_limits<double>::max();
          for (int j = m_or->indG[dep]; j < m_or->indG[dep + 1]; j++){

            if (m_or->nodeG[j] == second && m_or->wG[j] < w){
              w =  m_or->wG[j];
              edge_index = j;

            }
          }

       //   if (edge_index == -1) Rcpp::Rcout << "ok" << endl;
          m_result[edge_index] += dem[i];
          dep = second;
          pos += 1;
        }
      }


    }


    std::fill(Distances.begin(),Distances.end(),std::numeric_limits<double>::max());
    std::fill(parents.begin(),parents.end(), -1);
    std::fill(parents2.begin(),parents2.end(), -1);

  }
}

void aonGraphC::phastr(std::size_t begin, std::size_t end){
  DVec Distances(m_gr->nbnode, std::numeric_limits<double>::max());
  IVec parents(m_gr->nbnode, -1);
  IVec parents2(m_gr->nbnode, -1);

  for (std::size_t k=begin; k!=end;k++){



    int StartNode = m_od[k].first;

    Distances[StartNode] = 0.0;
    PQ Q;
    Q.push(std::make_pair(StartNode, 0.0));


    while (true) {

      if (Q.empty()){
        break;
      }

      if (!Q.empty()){

        int v=Q.top().first;
        double w=Q.top().second;
        Q.pop();

        if (w <= Distances[v]) {
          if (Stall_par(v, Distances, m_gr->nodeG, m_gr->wG, m_gr->indG)) continue;

          for (int i = m_gr->indGr[v]; i< m_gr->indGr[v+1]; i++){
            int v2 = m_gr->nodeGr[i];
            double w2 = m_gr->wGr[i];


            if (Distances[v] + w2 < Distances[v2]) {
              Distances[v2] = Distances[v] + w2;
              parents[v2] = v;
              Q.push(std::make_pair(v2, Distances[v2]));

            }
          }
        }

      }
    }

    //Backward

    for (int i=0; i < (m_gr->indG.size()-1); i++){

      for (int j=m_gr->indG[i]; j < m_gr->indG[i+1]; j++){


        if (Distances[m_gr->nodeG[j]]+m_gr->wG[j] < Distances[i]) {
          Distances[i] = Distances[m_gr->nodeG[j]]+m_gr->wG[j];
          parents2[i] = m_gr->nodeG[j];
        }

      }
    }

    IVec arr_vec = m_od[k].second.first;
    DVec dem = m_od[k].second.second;

    for (unsigned int i=0; i < arr_vec.size();i++){
      int end = arr_vec[i];

      IVec result2;
      int last_node = -1;
      for (auto p = end; p != -1; p = parents2[p]){
        result2.push_back(p);
        last_node = p;
      }


      reverse(result2.begin(),result2.end());

      if (last_node != -1){
        for (auto p = parents[last_node]; p != -1; p = parents[p]){
          result2.insert(result2.begin(), p);
          last_node = p;
        }
      }

      if (result2.size() > 1){
        //    m_gr->unpack(result2);
      }

      // load AON
      if (result2.size() > 1){

        int pos = 0;
        int dep = result2[pos];
        int second = result2[pos + 1];
        int arr = result2.back();


        while (dep != arr){

          second = result2[pos + 1];
          int edge_index = -1;
          double w = std::numeric_limits<double>::max();
          for (int j = m_or->indG[second]; j < m_or->indG[second + 1]; j++){

            if (m_or->nodeG[j] == dep && m_or->wG[j] < w){
              w =  m_or->wG[j];
              edge_index = j;

            }
          }

          //   if (edge_index == -1) Rcpp::Rcout << "ok" << endl;
          m_result[edge_index] += dem[i];
          dep = second;
          pos += 1;
        }
      }


    }


    std::fill(Distances.begin(),Distances.end(),std::numeric_limits<double>::max());
    std::fill(parents.begin(),parents.end(), -1);
    std::fill(parents2.begin(),parents2.end(), -1);

  }
}


unpackC::unpackC(CGraph* gr, Graph* original, Graph* cg, DVec input, bool convert):
  m_gr(gr), m_or(original), m_cg(cg), m_input(input), phast_rename(convert){
  m_result.resize(m_or->nbedge, 0.0);
  if (phast_rename){
    node_dict.resize(m_gr->nbnode);
    for (int i = 0; i < node_dict.size(); i++){
      node_dict[m_gr->nbnode - m_gr->rank[i]] = i;
    }
  }
}

unpackC::unpackC(unpackC &pn,  RcppParallel::Split):
  m_gr(pn.m_gr), m_or(pn.m_or), m_cg(pn.m_cg), m_input(pn.m_input), phast_rename(pn.phast_rename),
  m_result(pn.m_result), node_dict(pn.node_dict)
{
}

void unpackC::join(unpackC &pb){
  for (int i = 0; i < m_result.size(); i++){
    m_result[i] += pb.m_result[i];
  }
}

void unpackC::operator()(std::size_t begin, std::size_t end){

  for (size_t k = begin; k != end; k++){

    for (int i = m_cg->indG[k]; i < m_cg->indG[k+1]; i++){

      if (m_input[i] == 0) continue;
      IVec result2(2);
      result2[0] = k;
      result2[1] = m_cg->nodeG[i];

      m_gr->unpack(result2);

      int pos = 0;
      int dep = result2[pos];
      int second = result2[pos + 1];
      int arr = result2.back();

      if (phast_rename){
        dep = node_dict[dep];
        second = node_dict[second];
        arr = node_dict[arr];
      }

      while (dep != arr){

        second = result2[pos + 1];
        if (phast_rename){
          second = node_dict[second];
        }
        int edge_index = -1;
        double w = std::numeric_limits<double>::max();
        for (int j = m_or->indG[dep]; j < m_or->indG[dep + 1]; j++){

          if (m_or->nodeG[j] == second && m_or->wG[j] < w){
            w =  m_or->wG[j];
            edge_index = j;

          }
        }

        if (edge_index == -1) Rcpp::Rcout << dep << "->"<<second << endl;
     //   if (edge_index != -1) Rcpp::Rcout << "ok2" << endl;
        m_result[edge_index] += m_input[i];
        dep = second;
        pos += 1;
      }

    }

  }
}
