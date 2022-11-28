// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]

#include <vector>
#include <Rcpp.h>
#include <RcppParallel.h>
#include "graph.h"
#include "path_pair.h"
#include "stall.h"
using namespace RcppParallel;

// constructor
pathPair::pathPair(Graph *gr, IVec dep, IVec arr, IVec keep, double lim, int algo) :
  m_gr(gr), m_dep(dep), m_arr(arr), m_keep(keep), m_lim(lim), algorithm(algo)

{

  m_result.resize(m_dep.size());
}

// switch between algorithms
void pathPair::operator()(std::size_t begin, std::size_t end){
  if (algorithm == 0) dijkstra_early_stop(begin, end);
  if (algorithm == 1) bidir(begin, end);
  if (algorithm == 2) astar(begin, end);
  if (algorithm == 3) nba(begin, end);
  if (algorithm == 4) iso(begin, end);
  if (algorithm == 5) detour(begin, end);
}


// dijkstra early stop
void pathPair::dijkstra_early_stop(std::size_t begin, std::size_t end){

  std::vector<double> Distances(m_gr->nbnode, std::numeric_limits<double>::max());
  std::vector<int> parents(m_gr->nbnode, -1);

  for (std::size_t k = begin; k != end; k++){

    int StartNode=m_dep[k];

    Distances[StartNode] = 0.0;

    PQ Q;
    Q.push(std::make_pair(StartNode, 0.0));

    while (!Q.empty()) {
      int v = Q.top().first;
      double w = Q.top().second;
      Q.pop();

      if (w <= Distances[v]) {

        for (int i = m_gr->indG[v]; i<  m_gr->indG[v+1]; i++){
          int v2 = m_gr->nodeG[i];
          double w2 = m_gr->wG[i];

          if (Distances[v] + w2 < Distances[v2]) {
            Distances[v2] = Distances[v] + w2;
            parents[v2] = v;

            Q.push(std::make_pair(v2, Distances[v2]));
          }
        }
      }
      if (v == m_arr[k]){
        break;
      }
    }

    int end = m_arr[k];
    SVec result2;

    for (auto p = parents[end]; p != -1; p = parents[p]){
      if (m_keep[p] == 1) result2.push_back(m_gr->dict[p]);
    }

    if (result2.size()>0){
      if (m_keep[end] == 1) result2.insert(result2.begin(), m_gr->dict[end]);
      reverse(result2.begin(), result2.end());
    }

    m_result[k] = result2;


    //Reinitialize vectors
    std::fill(Distances.begin(),Distances.end(),std::numeric_limits<double>::max());
    std::fill(parents.begin(),parents.end(),-1);


  }

}


// bidirectional dijkstra
void pathPair::bidir(std::size_t begin, std::size_t end){

  std::vector<double> Distances(m_gr->nbnode, std::numeric_limits<double>::max());
  std::vector<double> Distances2(m_gr->nbnode, std::numeric_limits<double>::max());
  std::vector <int>  parents(m_gr->nbnode, -1);
  std::vector <int>  parents2(m_gr->nbnode, -1);
  std::vector <int> Visited(m_gr->nbnode,0);
  std::vector <int> Visiting(m_gr->nbnode,0);
  std::vector <int> Visited2(m_gr->nbnode,0);
  std::vector <int> Visiting2(m_gr->nbnode,0);

  for (std::size_t k=begin; k!=end;k++){

    int mid = -1;
    int StartNode = m_dep[k];
    int EndNode=m_arr[k];
    Distances[StartNode] = 0.0;
    Distances2[EndNode] = 0.0;

    PQ Q;
    PQ Qr;
    Q.push(std::make_pair(StartNode, 0.0));
    Qr.push(std::make_pair(EndNode, 0.0));
    Visiting[StartNode]=1;
    Visiting2[EndNode]=1;

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
          for (int i = m_gr->indG[v]; i< m_gr->indG[v+1]; i++){
            int v2 = m_gr->nodeG[i];
            double w2 = m_gr->wG[i];


            if (Distances[v] + w2 < Distances[v2]) {
              Distances[v2] = Distances[v] + w2;
              parents[v2] = v;

              Q.push(std::make_pair(v2, Distances[v2]));
              Visiting[v2]=1;

            }
          }
        }
        if ((Visited2[v]==1 || Visiting2[v]==1)  && (Distances[v]+Distances2[v]) < mu){

          mu=Distances[v]+Distances2[v];
          mid = v;

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
              parents2[vv2] = vv;

              Qr.push(std::make_pair(vv2, Distances2[vv2]));
              Visiting2[vv2]=1;
            }
          }
        }
        if ((Visited[vv]==1 || Visiting[vv]==1) && (Distances[vv]+Distances2[vv]) < mu){

          mu=Distances[vv]+Distances2[vv];
          mid = vv;
        }
      }

    }

    SVec result2;
    if (mid != -1){
      for (auto p = parents2[mid]; p != -1; p = parents2[p]){
        if (m_keep[p] == 1) result2.insert(result2.begin(),m_gr->dict[p]);
      }

      if (Distances[mid]!=numeric_limits<double>::max() || Distances2[mid]!=numeric_limits<double>::max()){
        if (m_keep[mid] == 1) result2.push_back(m_gr->dict[mid]);
      }

      for (auto p = parents[mid]; p != -1; p = parents[p]){
        if (m_keep[p] == 1) result2.push_back(m_gr->dict[p]);
      }

      reverse(result2.begin(), result2.end());
    }


    m_result[k] = result2;

    std::fill(Distances.begin(), Distances.end(), std::numeric_limits<double>::max());
    std::fill(Distances2.begin(), Distances2.end(), std::numeric_limits<double>::max());
    std::fill(parents.begin(), parents.end(), -1);
    std::fill(parents2.begin(), parents2.end(), -1);
    std::fill(Visited.begin(), Visited.end(), 0);
    std::fill(Visited2.begin(), Visited2.end(), 0);
    std::fill(Visiting.begin(), Visiting.end(), 0);
    std::fill(Visiting2.begin(), Visiting2.end(), 0);

  }
}


// A*
void pathPair::astar(std::size_t begin, std::size_t end){

  std::vector<double> Distances(m_gr->nbnode, std::numeric_limits<double>::max());
  std::vector<double> Distances2(m_gr->nbnode, std::numeric_limits<double>::max());
  std::vector<int> parents(m_gr->nbnode, -1);
  std::vector <int> closedList(m_gr->nbnode,0);
  std::vector <int> openList(m_gr->nbnode,0);
  for (std::size_t k=begin; k!=end;k++){

    int StartNode=m_dep[k];
    int endNode=m_arr[k];
    double lata=m_gr->lat[endNode];
    double lona=m_gr->lon[endNode];

    Distances[StartNode] = 0.0;
    Distances2[StartNode] = sqrt(pow(m_gr->lat[StartNode]-lata,2)+pow(m_gr->lon[StartNode]-lona,2))/m_gr->k;
    PQ Q;
    Q.push(std::make_pair(StartNode,sqrt(pow(m_gr->lat[StartNode]-lata,2)+pow(m_gr->lon[StartNode]-lona,2))/m_gr->k));
    openList[StartNode]=1;

    while (!Q.empty()) {
      int v = Q.top().first;
      Q.pop();
      if (closedList[v]==1){
        continue;
      }
      openList[v]=0;
      closedList[v]=1;

      for (int i=m_gr->indG[v]; i< m_gr->indG[v+1]; i++){
        int v2 = m_gr->nodeG[i];
        double w2 = m_gr->wG[i];
        if (closedList[v2]==1) {
          continue;
        }

        double temp;
        temp = Distances[v] + w2;
        if (openList[v2]==0){

          Q.push(std::make_pair(v2,Distances2[v2]));
          openList[v2]=1;
        }


        else if (temp>=Distances[v2]){
          continue;
        }


        Distances[v2]=temp;
        Distances2[v2]=Distances[v2]+sqrt(pow(m_gr->lat[v2]-lata,2)+pow(m_gr->lon[v2]-lona,2))/m_gr->k;
        parents[v2] = v;
        Q.push(std::make_pair(v2,Distances2[v2]));
        openList[v2]=1;
      }


      if (v==endNode){
        break;
      }

    }

    int end = m_arr[k];
    SVec result2;

    for (auto p = parents[end]; p != -1; p = parents[p]){
      if (m_keep[p] == 1) result2.push_back(m_gr->dict[p]);
    }

    if (result2.size()>0){
      if (m_keep[end] == 1) result2.insert(result2.begin(), m_gr->dict[end]);
      reverse(result2.begin(), result2.end());
    }

    m_result[k] = result2;

    //Reinitialize vectors
    std::fill(Distances.begin(),Distances.end(),std::numeric_limits<double>::max());
    std::fill(Distances2.begin(),Distances2.end(),std::numeric_limits<double>::max());
    std::fill(parents.begin(),parents.end(),-1);
    std::fill(closedList.begin(),closedList.end(),0);
    std::fill(openList.begin(),openList.end(),0);


  }
}



// bidirectional A*
void pathPair::nba(std::size_t begin, std::size_t end){

  std::vector<double> Distances(m_gr->nbnode, std::numeric_limits<double>::max());
  std::vector<double> Distances2(m_gr->nbnode, std::numeric_limits<double>::max());
  std::vector<int> parents(m_gr->nbnode, -1);
  std::vector<int> parents2(m_gr->nbnode, -1);
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



    //double Pr=0.5*sqrt(pow(lat[StartNode]-lata,2)+pow(lon[StartNode]-lona,2))/k;
    double total1=sqrt(pow(lata2-lata,2)+pow(lona2-lona,2))/m_gr->k;
    double total2=sqrt(pow(lata2-lata,2)+pow(lona2-lona,2))/m_gr->k;
    int mid;
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
                parents[v2]=v;

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
                parents2[vv2] = vv;

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

    SVec result2;
    if (mid > -1) {
      for (auto p = parents2[mid]; p != -1; p = parents2[p]){

        if (m_keep[p]==1)  result2.insert(result2.begin(),m_gr->dict[p]);
      }

      if (Distances[mid]!=numeric_limits<double>::max() || Distances2[mid]!=numeric_limits<double>::max()){
        if (m_keep[mid]==1) result2.push_back(m_gr->dict[mid]);
      }

      for (auto p = parents[mid]; p != -1; p = parents[p]){
        if (m_keep[p]==1) result2.push_back(m_gr->dict[p]);
      }

      reverse(result2.begin(), result2.end());
    }



    m_result[k] = result2;


    //Reinitialize
    std::fill(Distances.begin(),Distances.end(),std::numeric_limits<double>::max());
    std::fill(Distances2.begin(),Distances2.end(),std::numeric_limits<double>::max());
    std::fill(parents.begin(),parents.end(),-1);
    std::fill(parents2.begin(),parents2.end(),-1);
    std::fill(Visited.begin(),Visited.end(),0);
    std::fill(Visited1check.begin(),Visited1check.end(),0);
    std::fill(Visited2check.begin(),Visited2check.end(),0);

  }

}

// isochrone (one limit)
void pathPair::iso(std::size_t begin, std::size_t end){

  //Boucle sur chaque trajet
  std::vector<double> Distances(m_gr->nbnode, std::numeric_limits<double>::max());

  for (std::size_t k = begin; k != end; k++){

    int StartNode = m_dep[k];
    Distances[StartNode] = 0.0;
    PQ Q;
    Q.push(std::make_pair(StartNode, 0.0));

    while (!Q.empty()) {
      int v = Q.top().first;
      double w = Q.top().second;
      Q.pop();

      if (w <= Distances[v]) {


        for (int i = m_gr->indG[v]; i< m_gr->indG[v+1]; i++) {

          int v2 =  m_gr->nodeG[i];
          double w2 = m_gr->wG[i];

          if (Distances[v] + w2 < Distances[v2]) {
            Distances[v2] = Distances[v] + w2;

            Q.push(make_pair(v2, Distances[v2]));
          }

        }

      }
      if (Distances[v] > m_lim){
        break;
      }
    }

    SVec result2;

    for (int i = 0; i != Distances.size(); i++){
      if (Distances[i] < m_lim){
        if (m_keep[i]==1) result2.push_back(m_gr->dict[i]);
      }
    }
    m_result[k] = result2;


    //Reinitialize
    std::fill(Distances.begin(),Distances.end(),std::numeric_limits<double>::max());

  }
}


// detour

void pathPair::detour(std::size_t begin, std::size_t end){

  DVec distances(m_gr->nbnode, numeric_limits<double>::max());
  DVec distances2(m_gr->nbnode, numeric_limits<double>::max());
  IVec visited(m_gr->nbnode, 0);

  for (std::size_t k = begin; k != end; k++){


    int start = m_dep[k];
    int end = m_arr[k];

    distances[start] = 0.0;
    distances2[end] = 0.0;

    PQ Q;
    PQ Qr;
    Q.push(make_pair(start, 0.0));
    Qr.push(make_pair(end, 0.0));
    visited[start] += 1;
    visited[end] += 1;

    double mu=numeric_limits<double>::max();
    double def_mu=numeric_limits<double>::max();

    while (!Q.empty() && !Qr.empty()) {
      if (Q.top().second+Qr.top().second >= mu){
        def_mu=mu;
      }

      if (Q.top().second > def_mu + m_lim && Qr.top().second > def_mu + m_lim){
        break;
      }

      if (!Q.empty()){
        int v=Q.top().first;
        int w=Q.top().second;
        Q.pop();

        if (w <= distances[v]) {

          for (int i= m_gr->indG[v]; i< m_gr->indG[v+1]; i++){
            int v2 = m_gr->nodeG[i];
            double w2 = m_gr->wG[i];


            if (distances[v] + w2 < distances[v2]) {
              distances[v2] = distances[v] + w2;


              Q.push(make_pair(v2, distances[v2]));
              visited[v2]+=1;

            }
          }
        }

        if ((visited[v] > 1)  && (distances[v]+distances2[v]) < mu){

          mu=distances[v]+distances2[v];

        }
      }

      if (!Qr.empty()){
        int vv=Qr.top().first;
        int ww=Qr.top().second;
        Qr.pop();


        visited[vv]+=1;


        if (ww <= distances2[vv]) {
          for (int i = m_gr->indGr[vv]; i< m_gr->indGr[vv+1]; i++){
            int vv2 = m_gr->nodeGr[i];
            double ww2 = m_gr->wGr[i];


            if (distances2[vv] + ww2 < distances2[vv2]) {
              distances2[vv2] = distances2[vv] + ww2;


              Qr.push(make_pair(vv2, distances2[vv2]));
              visited[vv2]+=1;
            }
          }
        }


        if ((visited[vv]> 1) && (distances[vv]+distances2[vv]) < mu){

          mu=distances[vv]+distances2[vv];

        }
      }

    }


    if (mu >= numeric_limits<double>::max()){
      continue;
    }
    else {
      SVec result2;
      for (int i=0; i < distances.size(); ++i){
        if (distances[i]+distances2[i] < def_mu + m_lim){
          if (m_keep[i] == 1) result2.push_back(m_gr->dict[i]);
        }
      }
      m_result[k] = result2;
    }

    //Reinitialize
    fill(distances.begin(),distances.end(),numeric_limits<double>::max());
    fill(distances2.begin(),distances2.end(),numeric_limits<double>::max());
    fill(visited.begin(),visited.end(),0);
  }
}


////////////////////////////////////////////////// contracted graph
// constructor

pathPairC::pathPairC(CGraph *gr, IVec dep, IVec arr, IVec keep, int algo) :
  m_gr(gr), m_dep(dep), m_arr(arr), m_keep(keep), algorithm(algo)
{
  m_result.resize(m_dep.size());
}


void pathPairC::operator()(std::size_t begin, std::size_t end){
  if (algorithm == 0) bidirmod(begin, end);

}

void pathPairC::bidirmod(std::size_t begin, std::size_t end){
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

        m_gr->unpack(result2);

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


    SVec result3;
    for (int i = 0; i < result2.size(); i++) if (m_keep[result2[i]]==1) result3.push_back(m_gr->dict[result2[i]]);

    m_result[k]=result3;

  }

}

