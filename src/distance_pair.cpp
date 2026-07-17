// [[Rcpp::depends(RcppParallel)]]


#include <vector>
#include <Rcpp.h>
#include <RcppParallel.h>
#include "graph.h"
#include "distance_pair.h"
#include "stall.h"
using namespace RcppParallel;
//////////////////////////////////////////////////////////////// Normal graph
// constructor
distancePair::distancePair(Graph *gr, IVec dep, IVec arr, int algo) :
  m_gr(gr), m_dep(dep), m_arr(arr), algorithm(algo)
{
  add = false;
  if (m_gr->add.size() > 0) add = true;
  m_result.resize(m_dep.size(), 0.0);
}

// switch between algorigthms
void distancePair::operator()(std::size_t begin, std::size_t end){
  if (algorithm == 0) dijkstra_early_stop(begin, end);
  if (algorithm == 1) bidir(begin, end);
  if (algorithm == 2) astar(begin, end);
  if (algorithm == 3) nba(begin, end);
}


// dijkstra early stop
void distancePair::dijkstra_early_stop(std::size_t begin, std::size_t end){

  std::vector<double> Distances(m_gr->nbnode, std::numeric_limits<double>::max());
  DVec distances2;
  if (add) distances2.resize(m_gr->nbnode, numeric_limits<double>::max());

  for (std::size_t k = begin; k != end; k++){

    int StartNode=m_dep[k];

    Distances[StartNode] = 0.0;
    if (add) distances2[StartNode] = 0.0;

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
            if (add) distances2[v2] = distances2[v] + m_gr->add[i];

            Q.push(std::make_pair(v2, Distances[v2]));
          }
        }
      }
      if (v == m_arr[k]){
        break;
      }
    }

    if (add){
      m_result[k] = distances2[m_arr[k]];
    } else{
      m_result[k] = Distances[m_arr[k]];
    }



    //Reinitialize vectors
    std::fill(Distances.begin(),Distances.end(),std::numeric_limits<double>::max());
    if (add) std::fill(distances2.begin(),distances2.end(), std::numeric_limits<double>::max());

  }

}

// bidirectional dijkstra
void distancePair::bidir(std::size_t begin, std::size_t end){

  std::vector<double> Distances(m_gr->nbnode, std::numeric_limits<double>::max());
  std::vector<double> Distances2(m_gr->nbnode, std::numeric_limits<double>::max());
  DVec distadd;
  DVec distadd2;
  if (add){
    distadd.resize(m_gr->nbnode, std::numeric_limits<double>::max());
    distadd2.resize(m_gr->nbnode, std::numeric_limits<double>::max());
  }
  std::vector <int> Visited(m_gr->nbnode,0);
  std::vector <int> Visiting(m_gr->nbnode,0);
  std::vector <int> Visited2(m_gr->nbnode,0);
  std::vector <int> Visiting2(m_gr->nbnode,0);

  for (std::size_t k=begin; k!=end;k++){

    int StartNode = m_dep[k];
    int EndNode=m_arr[k];
    Distances[StartNode] = 0.0;
    Distances2[EndNode] = 0.0;
    if (add){
      distadd[StartNode] = 0.0;
      distadd2[EndNode] = 0.0;
    }

    PQ Q;
    PQ Qr;
    Q.push(std::make_pair(StartNode, 0.0));
    Qr.push(std::make_pair(EndNode, 0.0));
    Visiting[StartNode]=1;
    Visiting2[EndNode]=1;

    double mu=std::numeric_limits<double>::max();
    double mu2 = std::numeric_limits<double>::max(); // additional weight

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
              if (add) distadd[v2] = distadd[v] + m_gr->add[i];

              Q.push(std::make_pair(v2, Distances[v2]));
              Visiting[v2]=1;

            }
          }
        }
        if ((Visited2[v]==1 || Visiting2[v]==1)  && (Distances[v]+Distances2[v]) < mu){

          mu=Distances[v]+Distances2[v];
          if (add)  mu2 = distadd[v]+distadd2[v];

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
              if (add) distadd2[vv2] = distadd2[vv] + m_gr->addr[i];

              Qr.push(std::make_pair(vv2, Distances2[vv2]));
              Visiting2[vv2]=1;
            }
          }
        }
        if ((Visited[vv]==1 || Visiting[vv]==1) && (Distances[vv]+Distances2[vv]) < mu){

          mu=Distances[vv]+Distances2[vv];
          if (add)  mu2 = distadd[vv]+distadd2[vv];

        }
      }

    }

    if (add){
      m_result[k]=mu2;
    } else{
      m_result[k]=mu;
    }

    std::fill(Distances.begin(), Distances.end(), std::numeric_limits<double>::max());
    std::fill(Distances2.begin(), Distances2.end(), std::numeric_limits<double>::max());
    if (add){
      std::fill(distadd.begin(), distadd.end(), std::numeric_limits<double>::max());
      std::fill(distadd2.begin(), distadd2.end(), std::numeric_limits<double>::max());
    }
    std::fill(Visited.begin(), Visited.end(), 0);
    std::fill(Visited2.begin(), Visited2.end(), 0);
    std::fill(Visiting.begin(), Visiting.end(), 0);
    std::fill(Visiting2.begin(), Visiting2.end(), 0);

  }
}

// A*
void distancePair::astar(std::size_t begin, std::size_t end){

  std::vector<double> Distances(m_gr->nbnode, std::numeric_limits<double>::max());
  DVec distadd;
  if (add){
    distadd.resize(m_gr->nbnode, std::numeric_limits<double>::max());
  }
  std::vector<double> Distances2(m_gr->nbnode, std::numeric_limits<double>::max());
  std::vector <int> closedList(m_gr->nbnode,0);
  std::vector <int> openList(m_gr->nbnode,0);
  for (std::size_t k=begin; k!=end;k++){

    int StartNode=m_dep[k];
    int endNode=m_arr[k];
    double lata=m_gr->lat[endNode];
    double lona=m_gr->lon[endNode];

    Distances[StartNode] = 0.0;
    if (add) distadd[StartNode] = 0.0;

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
        if (add) distadd[v2] = distadd[v] + m_gr->add[i];

        Distances2[v2]=Distances[v2]+sqrt(pow(m_gr->lat[v2]-lata,2)+pow(m_gr->lon[v2]-lona,2))/m_gr->k;
        Q.push(std::make_pair(v2,Distances2[v2]));
        openList[v2]=1;
      }


      if (v==endNode){
        break;
      }

    }

    if (add){
      m_result[k]= distadd[endNode];
    } else{
      m_result[k]= Distances[endNode];
    }


    //Reinitialize vectors
    std::fill(Distances.begin(),Distances.end(),std::numeric_limits<double>::max());
    if (add) std::fill(distadd.begin(),distadd.end(),std::numeric_limits<double>::max());
    std::fill(Distances2.begin(),Distances2.end(),std::numeric_limits<double>::max());
    std::fill(closedList.begin(),closedList.end(),0);
    std::fill(openList.begin(),openList.end(),0);


  }
}

// bidirectional A*
void distancePair::nba(std::size_t begin, std::size_t end){

  std::vector<double> Distances(m_gr->nbnode, std::numeric_limits<double>::max());
  std::vector<double> Distances2(m_gr->nbnode, std::numeric_limits<double>::max());
  DVec distadd;
  DVec distadd2;
  if (add){
    distadd.resize(m_gr->nbnode, std::numeric_limits<double>::max());
    distadd2.resize(m_gr->nbnode, std::numeric_limits<double>::max());
  }
  std::vector <int> Visited(m_gr->nbnode,0);
  std::vector <int> Visited1check(m_gr->nbnode,0);
  std::vector <int> Visited2check(m_gr->nbnode,0);

  for (std::size_t k=begin; k!=end;k++){

    int StartNode=m_dep[k];
    int EndNode=m_arr[k];

    if (StartNode==EndNode){
      m_result[k]=0;
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
    if (add){

      distadd[StartNode] = 0.0;
      distadd2[EndNode] = 0.0;
    }
    PQ Q;
    PQ Qr;
    Q.push(std::make_pair(StartNode, sqrt(pow(m_gr->lat[StartNode]-lata,2)+pow(m_gr->lon[StartNode]-lona,2))/m_gr->k));
    Qr.push(std::make_pair(EndNode, sqrt(pow(m_gr->lat[EndNode]-lata2,2)+pow(m_gr->lon[EndNode]-lona2,2))/m_gr->k));



    double total1=sqrt(pow(lata2-lata,2)+pow(lona2-lona,2))/m_gr->k;
    double total2=sqrt(pow(lata2-lata,2)+pow(lona2-lona,2))/m_gr->k;
    int mid;
    double mu=std::numeric_limits<double>::max();
    double mu2=std::numeric_limits<double>::max();

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
                if (add) distadd[v2] = distadd[v] + m_gr->add[i];
                Visited1check[v2]=1;

                Q.push(std::make_pair(v2, tentative + sqrt(pow(m_gr->lat[v2]-lata,2)+pow(m_gr->lon[v2]-lona,2))/m_gr->k));

                if (Visited2check[v2]==1){
                  double temp=tentative + Distances2[v2];
                  if (mu > temp){
                    mu=temp;
                    mid=v2;
                    if (add) mu2 = distadd[v] + m_gr->add[i] + distadd2[v2];

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
                if (add) distadd2[vv2] = distadd2[vv] + m_gr->addr[i];
                Visited2check[vv2]=1;

                Qr.push(std::make_pair(vv2, tentative + sqrt(pow(m_gr->lat[vv2]-lata2,2)+pow(m_gr->lon[vv2]-lona2,2))/m_gr->k));

                if (Visited1check[vv2]==1){
                  double temp=tentative + Distances[vv2];
                  if (mu > temp){
                    mu=temp;
                    mid=vv2;
                    if (add) mu2 = distadd2[vv] + m_gr->addr[i] + distadd[vv2];

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


    if (add){
      m_result[k]=mu2;
    } else{
      m_result[k]=mu;
    }

    //Reinitialize
    std::fill(Distances.begin(),Distances.end(),std::numeric_limits<double>::max());
    std::fill(Distances2.begin(),Distances2.end(),std::numeric_limits<double>::max());
    if (add){
      std::fill(distadd.begin(),distadd.end(),std::numeric_limits<double>::max());
      std::fill(distadd2.begin(),distadd2.end(),std::numeric_limits<double>::max());
    }
    std::fill(Visited.begin(),Visited.end(),0);
    std::fill(Visited1check.begin(),Visited1check.end(),0);
    std::fill(Visited2check.begin(),Visited2check.end(),0);

  }

}


//////////////////////////////////////////////////////////// Normal graph, N aux columns

// constructor
distancePairAdd::distancePairAdd(Graph *gr, IVec dep, IVec arr, int algo) :
  m_gr(gr), m_dep(dep), m_arr(arr), algorithm(algo)
{
  n_add = m_gr->n_add;
  m_result.assign((size_t)m_dep.size() * n_add, 0.0);
}

void distancePairAdd::operator()(std::size_t begin, std::size_t end){
  if (algorithm == 0) dijkstra_early_stop(begin, end);
  if (algorithm == 1) bidir(begin, end);
  if (algorithm == 2) astar(begin, end);
  if (algorithm == 3) nba(begin, end);
}

// dijkstra early stop, N aux columns summed inline (same structure as
// distancePair::dijkstra_early_stop, generalized to n_add running sums, with
// touched-array resets instead of full fills)
void distancePairAdd::dijkstra_early_stop(std::size_t begin, std::size_t end){

  int n_pairs = (int)m_dep.size();
  int total = (int)m_gr->nodeG.size();

  std::vector<double> Distances(m_gr->nbnode, std::numeric_limits<double>::max());
  DVec distadd((size_t)m_gr->nbnode * n_add, std::numeric_limits<double>::max());
  IVec touched;

  for (std::size_t k = begin; k != end; k++){

    int StartNode=m_dep[k];
    int EndNode=m_arr[k];

    Distances[StartNode] = 0.0;
    for (int col = 0; col < n_add; col++){
      distadd[(size_t)col * m_gr->nbnode + StartNode] = 0.0;
    }
    touched.push_back(StartNode);

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
            for (int col = 0; col < n_add; col++){
              distadd[(size_t)col * m_gr->nbnode + v2] = distadd[(size_t)col * m_gr->nbnode + v] + m_gr->add[(size_t)col * total + i];
            }

            Q.push(std::make_pair(v2, Distances[v2]));
            touched.push_back(v2);
          }
        }
      }
      if (v == EndNode){
        break;
      }
    }

    for (int col = 0; col < n_add; col++){
      m_result[(size_t)col * n_pairs + k] = distadd[(size_t)col * m_gr->nbnode + EndNode];
    }

    for (int i = 0; i < (int)touched.size(); i++){
      int node = touched[i];
      Distances[node] = std::numeric_limits<double>::max();
      for (int col = 0; col < n_add; col++){
        distadd[(size_t)col * m_gr->nbnode + node] = std::numeric_limits<double>::max();
      }
    }
    touched.clear();

  }

}

// bidirectional dijkstra, N aux columns summed inline (same structure as
// distancePair::bidir, generalized from one running sum to n_add of them)
void distancePairAdd::bidir(std::size_t begin, std::size_t end){

  int n_pairs = (int)m_dep.size();
  int total = (int)m_gr->nodeG.size(); // stride of m_gr->add
  int totalr = (int)m_gr->nodeGr.size(); // stride of m_gr->addr

  std::vector<double> Distances(m_gr->nbnode, std::numeric_limits<double>::max());
  std::vector<double> Distances2(m_gr->nbnode, std::numeric_limits<double>::max());
  DVec distadd((size_t)m_gr->nbnode * n_add, std::numeric_limits<double>::max());
  DVec distadd2((size_t)m_gr->nbnode * n_add, std::numeric_limits<double>::max());
  DVec mu2(n_add, std::numeric_limits<double>::max());

  std::vector <int> Visited(m_gr->nbnode,0);
  std::vector <int> Visiting(m_gr->nbnode,0);
  std::vector <int> Visited2(m_gr->nbnode,0);
  std::vector <int> Visiting2(m_gr->nbnode,0);
  IVec touched;

  for (std::size_t k=begin; k!=end;k++){

    int StartNode = m_dep[k];
    int EndNode=m_arr[k];
    Distances[StartNode] = 0.0;
    Distances2[EndNode] = 0.0;
    for (int col = 0; col < n_add; col++){
      distadd[(size_t)col * m_gr->nbnode + StartNode] = 0.0;
      distadd2[(size_t)col * m_gr->nbnode + EndNode] = 0.0;
    }
    touched.push_back(StartNode);
    touched.push_back(EndNode);

    PQ Q;
    PQ Qr;
    Q.push(std::make_pair(StartNode, 0.0));
    Qr.push(std::make_pair(EndNode, 0.0));
    Visiting[StartNode]=1;
    Visiting2[EndNode]=1;

    double mu=std::numeric_limits<double>::max();
    std::fill(mu2.begin(), mu2.end(), std::numeric_limits<double>::max());

    while (!Q.empty() && !Qr.empty()) {
      if (Q.top().second+Qr.top().second >= mu){
        break;
      }


      if (!Q.empty()){
        int v=Q.top().first;
        double w=Q.top().second;
        Q.pop();


        Visited[v]=1;


        if (w <= Distances[v]) {
          for (int i = m_gr->indG[v]; i< m_gr->indG[v+1]; i++){
            int v2 = m_gr->nodeG[i];
            double w2 = m_gr->wG[i];


            if (Distances[v] + w2 < Distances[v2]) {
              Distances[v2] = Distances[v] + w2;
              for (int col = 0; col < n_add; col++){
                distadd[(size_t)col * m_gr->nbnode + v2] = distadd[(size_t)col * m_gr->nbnode + v] + m_gr->add[(size_t)col * total + i];
              }

              Q.push(std::make_pair(v2, Distances[v2]));
              Visiting[v2]=1;
              touched.push_back(v2);

            }
          }
        }
        if ((Visited2[v]==1 || Visiting2[v]==1)  && (Distances[v]+Distances2[v]) < mu){

          mu=Distances[v]+Distances2[v];
          for (int col = 0; col < n_add; col++){
            mu2[col] = distadd[(size_t)col * m_gr->nbnode + v] + distadd2[(size_t)col * m_gr->nbnode + v];
          }

        }
      }

      if (!Qr.empty()){
        int vv=Qr.top().first;
        double ww=Qr.top().second;
        Qr.pop();


        Visited2[vv]=1;

        if (ww <= Distances2[vv]) {
          for (int i=m_gr->indGr[vv]; i< m_gr->indGr[vv+1]; i++){
            int vv2 = m_gr->nodeGr[i];
            double ww2 = m_gr->wGr[i];


            if (Distances2[vv] + ww2 < Distances2[vv2]) {
              Distances2[vv2] = Distances2[vv] + ww2;
              for (int col = 0; col < n_add; col++){
                distadd2[(size_t)col * m_gr->nbnode + vv2] = distadd2[(size_t)col * m_gr->nbnode + vv] + m_gr->addr[(size_t)col * totalr + i];
              }

              Qr.push(std::make_pair(vv2, Distances2[vv2]));
              Visiting2[vv2]=1;
              touched.push_back(vv2);
            }
          }
        }
        if ((Visited[vv]==1 || Visiting[vv]==1) && (Distances[vv]+Distances2[vv]) < mu){

          mu=Distances[vv]+Distances2[vv];
          for (int col = 0; col < n_add; col++){
            mu2[col] = distadd[(size_t)col * m_gr->nbnode + vv] + distadd2[(size_t)col * m_gr->nbnode + vv];
          }

        }
      }

    }

    for (int col = 0; col < n_add; col++){
      m_result[(size_t)col * n_pairs + k] = mu2[col];
    }

    for (int i = 0; i < touched.size(); i++){
      int node = touched[i];
      Distances[node] = std::numeric_limits<double>::max();
      Distances2[node] = std::numeric_limits<double>::max();
      for (int col = 0; col < n_add; col++){
        distadd[(size_t)col * m_gr->nbnode + node] = std::numeric_limits<double>::max();
        distadd2[(size_t)col * m_gr->nbnode + node] = std::numeric_limits<double>::max();
      }
      Visited[node] = 0;
      Visited2[node] = 0;
      Visiting[node] = 0;
      Visiting2[node] = 0;
    }
    touched.clear();

  }
}

// A*, N aux columns summed inline (same structure as distancePair::astar,
// generalized to n_add running sums, with touched-array resets instead of
// full fills)
void distancePairAdd::astar(std::size_t begin, std::size_t end){

  int n_pairs = (int)m_dep.size();
  int total = (int)m_gr->nodeG.size();

  std::vector<double> Distances(m_gr->nbnode, std::numeric_limits<double>::max());
  DVec distadd((size_t)m_gr->nbnode * n_add, std::numeric_limits<double>::max());
  std::vector<double> Distances2(m_gr->nbnode, std::numeric_limits<double>::max());
  std::vector <int> closedList(m_gr->nbnode,0);
  std::vector <int> openList(m_gr->nbnode,0);
  IVec touched;

  for (std::size_t k=begin; k!=end;k++){

    int StartNode=m_dep[k];
    int endNode=m_arr[k];
    double lata=m_gr->lat[endNode];
    double lona=m_gr->lon[endNode];

    Distances[StartNode] = 0.0;
    for (int col = 0; col < n_add; col++){
      distadd[(size_t)col * m_gr->nbnode + StartNode] = 0.0;
    }

    Distances2[StartNode] = sqrt(pow(m_gr->lat[StartNode]-lata,2)+pow(m_gr->lon[StartNode]-lona,2))/m_gr->k;
    PQ Q;
    Q.push(std::make_pair(StartNode,Distances2[StartNode]));
    openList[StartNode]=1;
    touched.push_back(StartNode);

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
          touched.push_back(v2);
        }


        else if (temp>=Distances[v2]){
          continue;
        }


        Distances[v2]=temp;
        for (int col = 0; col < n_add; col++){
          distadd[(size_t)col * m_gr->nbnode + v2] = distadd[(size_t)col * m_gr->nbnode + v] + m_gr->add[(size_t)col * total + i];
        }

        Distances2[v2]=Distances[v2]+sqrt(pow(m_gr->lat[v2]-lata,2)+pow(m_gr->lon[v2]-lona,2))/m_gr->k;
        Q.push(std::make_pair(v2,Distances2[v2]));
        openList[v2]=1;
      }


      if (v==endNode){
        break;
      }

    }

    for (int col = 0; col < n_add; col++){
      m_result[(size_t)col * n_pairs + k] = distadd[(size_t)col * m_gr->nbnode + endNode];
    }

    for (int i = 0; i < (int)touched.size(); i++){
      int node = touched[i];
      Distances[node] = std::numeric_limits<double>::max();
      for (int col = 0; col < n_add; col++){
        distadd[(size_t)col * m_gr->nbnode + node] = std::numeric_limits<double>::max();
      }
      Distances2[node] = std::numeric_limits<double>::max();
      closedList[node] = 0;
      openList[node] = 0;
    }
    touched.clear();

  }
}

// bidirectional A*, N aux columns summed inline (same structure as
// distancePair::nba, generalized to n_add running sums, with touched-array
// resets instead of full fills)
void distancePairAdd::nba(std::size_t begin, std::size_t end){

  int n_pairs = (int)m_dep.size();
  int total = (int)m_gr->nodeG.size();
  int totalr = (int)m_gr->nodeGr.size();

  std::vector<double> Distances(m_gr->nbnode, std::numeric_limits<double>::max());
  std::vector<double> Distances2(m_gr->nbnode, std::numeric_limits<double>::max());
  DVec distadd((size_t)m_gr->nbnode * n_add, std::numeric_limits<double>::max());
  DVec distadd2((size_t)m_gr->nbnode * n_add, std::numeric_limits<double>::max());
  DVec mu2(n_add, std::numeric_limits<double>::max());
  std::vector <int> Visited(m_gr->nbnode,0);
  std::vector <int> Visited1check(m_gr->nbnode,0);
  std::vector <int> Visited2check(m_gr->nbnode,0);
  IVec touched;

  for (std::size_t k=begin; k!=end;k++){

    int StartNode=m_dep[k];
    int EndNode=m_arr[k];

    if (StartNode==EndNode){
      for (int col = 0; col < n_add; col++){
        m_result[(size_t)col * n_pairs + k] = 0;
      }
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
    for (int col = 0; col < n_add; col++){
      distadd[(size_t)col * m_gr->nbnode + StartNode] = 0.0;
      distadd2[(size_t)col * m_gr->nbnode + EndNode] = 0.0;
    }
    touched.push_back(StartNode);
    touched.push_back(EndNode);

    PQ Q;
    PQ Qr;
    Q.push(std::make_pair(StartNode, sqrt(pow(m_gr->lat[StartNode]-lata,2)+pow(m_gr->lon[StartNode]-lona,2))/m_gr->k));
    Qr.push(std::make_pair(EndNode, sqrt(pow(m_gr->lat[EndNode]-lata2,2)+pow(m_gr->lon[EndNode]-lona2,2))/m_gr->k));



    double total1=sqrt(pow(lata2-lata,2)+pow(lona2-lona,2))/m_gr->k;
    double total2=sqrt(pow(lata2-lata,2)+pow(lona2-lona,2))/m_gr->k;
    double mu=std::numeric_limits<double>::max();
    std::fill(mu2.begin(), mu2.end(), std::numeric_limits<double>::max());

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
                for (int col = 0; col < n_add; col++){
                  distadd[(size_t)col * m_gr->nbnode + v2] = distadd[(size_t)col * m_gr->nbnode + v] + m_gr->add[(size_t)col * total + i];
                }
                if (Visited1check[v2]==0) touched.push_back(v2);
                Visited1check[v2]=1;

                Q.push(std::make_pair(v2, tentative + sqrt(pow(m_gr->lat[v2]-lata,2)+pow(m_gr->lon[v2]-lona,2))/m_gr->k));

                if (Visited2check[v2]==1){
                  double temp=tentative + Distances2[v2];
                  if (mu > temp){
                    mu=temp;
                    for (int col = 0; col < n_add; col++){
                      mu2[col] = distadd[(size_t)col * m_gr->nbnode + v] + m_gr->add[(size_t)col * total + i] + distadd2[(size_t)col * m_gr->nbnode + v2];
                    }
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
                for (int col = 0; col < n_add; col++){
                  distadd2[(size_t)col * m_gr->nbnode + vv2] = distadd2[(size_t)col * m_gr->nbnode + vv] + m_gr->addr[(size_t)col * totalr + i];
                }
                if (Visited2check[vv2]==0) touched.push_back(vv2);
                Visited2check[vv2]=1;

                Qr.push(std::make_pair(vv2, tentative + sqrt(pow(m_gr->lat[vv2]-lata2,2)+pow(m_gr->lon[vv2]-lona2,2))/m_gr->k));

                if (Visited1check[vv2]==1){
                  double temp=tentative + Distances[vv2];
                  if (mu > temp){
                    mu=temp;
                    for (int col = 0; col < n_add; col++){
                      mu2[col] = distadd2[(size_t)col * m_gr->nbnode + vv] + m_gr->addr[(size_t)col * totalr + i] + distadd[(size_t)col * m_gr->nbnode + vv2];
                    }
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


    for (int col = 0; col < n_add; col++){
      m_result[(size_t)col * n_pairs + k] = mu2[col];
    }

    for (int i = 0; i < (int)touched.size(); i++){
      int node = touched[i];
      Distances[node] = std::numeric_limits<double>::max();
      Distances2[node] = std::numeric_limits<double>::max();
      for (int col = 0; col < n_add; col++){
        distadd[(size_t)col * m_gr->nbnode + node] = std::numeric_limits<double>::max();
        distadd2[(size_t)col * m_gr->nbnode + node] = std::numeric_limits<double>::max();
      }
      Visited[node] = 0;
      Visited1check[node] = 0;
      Visited2check[node] = 0;
    }
    touched.clear();

  }

}


//////////////////////////////////////////////////////////////// Contracted graph

// constructor
distancePairC::distancePairC(CGraph *gr, IVec dep, IVec arr, int algo) :
  m_gr(gr), m_dep(dep), m_arr(arr), algorithm(algo)
{
  add = false;
  if (m_gr->add.size() > 0) add = true;
  m_result.resize(m_dep.size(), 0.0);
}


// overload operator()
void distancePairC::operator()(std::size_t begin, std::size_t end){
  if (algorithm == 0) bidirmod(begin, end);
}

// modified bidirectional dijkstra
void distancePairC::bidirmod(std::size_t begin, std::size_t end){

  std::vector<double> Distances(m_gr->nbnode, std::numeric_limits<double>::max());
  std::vector<double> Distances2(m_gr->nbnode, std::numeric_limits<double>::max());
  DVec distadd;
  DVec distadd2;


  if (add){
    distadd.resize(m_gr->nbnode, std::numeric_limits<double>::max());
    distadd2.resize(m_gr->nbnode, std::numeric_limits<double>::max());
  }

  std::vector<int> Visited1(m_gr->nbnode, 0);
  std::vector<int> Visited2(m_gr->nbnode, 0);
  std::vector<int> Visited;

  double mu=std::numeric_limits<double>::max();
  double mu2=std::numeric_limits<double>::max();


  for (std::size_t k=begin; k!=end;k++){

    int StartNode = m_dep [k];
    int EndNode = m_arr [k];

    Distances[StartNode] = 0.0;
    Distances2[EndNode] = 0.0;

    if (add){
      distadd[StartNode] = 0.0;
      distadd2[EndNode] = 0.0;
    }
    PQ Q;
    PQ Qr;
    Q.push(std::make_pair(StartNode, 0.0));
    Qr.push(std::make_pair(EndNode, 0.0));


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


        Visited1[v] = 1;

        if ((Visited2[v]==1)  && (Distances[v]+Distances2[v]) < mu){

          mu=Distances[v]+Distances2[v];
          if (add)  mu2 = distadd[v]+distadd2[v];

        }

        if (w <= Distances[v]) {
          if (Stall_par(v, Distances, m_gr->nodeGr, m_gr->wGr, m_gr->indGr)) continue;
          for (int i=m_gr->indG[v]; i< m_gr->indG[v+1]; i++){
            int v2 = m_gr->nodeG[i];
            double w2 = m_gr->wG[i];

            if (Distances[v] + w2 < Distances[v2]) {
              Distances[v2] = Distances[v] + w2;
              if (add) distadd[v2] = distadd[v] + m_gr->add[i];
              Q.push(std::make_pair(v2, Distances[v2]));
              Visited1[v2]=1;

              Visited.push_back(v2);

            }
          }
        }

      }

      if (!Qr.empty()){
        int vv = Qr.top().first;
        double ww = Qr.top().second;
        Qr.pop();

        Visited.push_back(vv);


        Visited2[vv]=1;


        if ((Visited1[vv]== 1) && (Distances[vv]+Distances2[vv]) < mu){

          mu=Distances[vv]+Distances2[vv];
          if (add)  mu2 = distadd[vv]+distadd2[vv];

        }

        if (ww <= Distances2[vv]) {
          if (Stall_par(vv, Distances2, m_gr->nodeG, m_gr->wG, m_gr->indG)) continue;

          for (int i=m_gr->indGr[vv]; i< m_gr->indGr[vv+1]; i++){
            int vv2 = m_gr->nodeGr[i];
            double ww2 = m_gr->wGr[i];


            if (Distances2[vv] + ww2 < Distances2[vv2]) {
              Distances2[vv2] = Distances2[vv] + ww2;
              if (add) distadd2[vv2] = distadd2[vv] + m_gr->addr[i];
              Qr.push(std::make_pair(vv2, Distances2[vv2]));
              Visited2[vv2]=1;

              Visited.push_back(vv2);

            }
          }
        }

      }

    }


    if (add){
      m_result[k]=mu2;
    } else{

      m_result[k]=mu;
    }

    for (int i=0; i<Visited.size();i++) {
      Distances[Visited[i]] = std::numeric_limits<double>::max();

      if (add){
        distadd[Visited[i]] = std::numeric_limits<double>::max();
        distadd2[Visited[i]] = std::numeric_limits<double>::max();
      }
    }

    for (int i=0; i<Visited.size();i++) Distances2[Visited[i]] = std::numeric_limits<double>::max();
    for (int i=0; i<Visited.size();i++) Visited1[Visited[i]] = 0;
    for (int i=0; i<Visited.size();i++) Visited2[Visited[i]] = 0;

    Visited.clear();

    mu=std::numeric_limits<double>::max();
    mu2=std::numeric_limits<double>::max();
  }

}


//////////////////////////////////////////////////////////// Contracted graph, N aux columns

// constructor
distancePairCAdd::distancePairCAdd(CGraph *gr, IVec dep, IVec arr) :
  m_gr(gr), m_dep(dep), m_arr(arr)
{
  n_add = m_gr->n_add;
  m_result.assign((size_t)m_dep.size() * n_add, 0.0);
}

void distancePairCAdd::operator()(std::size_t begin, std::size_t end){
  bidirmod(begin, end);
}

// modified bidirectional dijkstra, N aux columns summed inline (same
// structure as distancePairC::bidirmod, generalized from one running sum to
// n_add of them; add/addr already hold per-shortcut sums precomputed by
// aggCAdd, so no unpacking happens here)
void distancePairCAdd::bidirmod(std::size_t begin, std::size_t end){

  int n_pairs = (int)m_dep.size();
  int total = (int)m_gr->nodeG.size(); // stride of m_gr->add
  int totalr = (int)m_gr->nodeGr.size(); // stride of m_gr->addr

  std::vector<double> Distances(m_gr->nbnode, std::numeric_limits<double>::max());
  std::vector<double> Distances2(m_gr->nbnode, std::numeric_limits<double>::max());
  DVec distadd((size_t)m_gr->nbnode * n_add, std::numeric_limits<double>::max());
  DVec distadd2((size_t)m_gr->nbnode * n_add, std::numeric_limits<double>::max());
  DVec mu2(n_add, std::numeric_limits<double>::max());

  std::vector<int> Visited1(m_gr->nbnode, 0);
  std::vector<int> Visited2(m_gr->nbnode, 0);
  std::vector<int> touched;

  double mu=std::numeric_limits<double>::max();

  for (std::size_t k=begin; k!=end;k++){

    int StartNode = m_dep[k];
    int EndNode = m_arr[k];

    Distances[StartNode] = 0.0;
    Distances2[EndNode] = 0.0;
    for (int col = 0; col < n_add; col++){
      distadd[(size_t)col * m_gr->nbnode + StartNode] = 0.0;
      distadd2[(size_t)col * m_gr->nbnode + EndNode] = 0.0;
    }
    touched.push_back(StartNode);
    touched.push_back(EndNode);

    PQ Q;
    PQ Qr;
    Q.push(std::make_pair(StartNode, 0.0));
    Qr.push(std::make_pair(EndNode, 0.0));


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

        touched.push_back(v);

        Visited1[v] = 1;

        if ((Visited2[v]==1)  && (Distances[v]+Distances2[v]) < mu){

          mu=Distances[v]+Distances2[v];
          for (int col = 0; col < n_add; col++){
            mu2[col] = distadd[(size_t)col * m_gr->nbnode + v] + distadd2[(size_t)col * m_gr->nbnode + v];
          }

        }

        if (w <= Distances[v]) {
          if (Stall_par(v, Distances, m_gr->nodeGr, m_gr->wGr, m_gr->indGr)) continue;
          for (int i=m_gr->indG[v]; i< m_gr->indG[v+1]; i++){
            int v2 = m_gr->nodeG[i];
            double w2 = m_gr->wG[i];

            if (Distances[v] + w2 < Distances[v2]) {
              Distances[v2] = Distances[v] + w2;
              for (int col = 0; col < n_add; col++){
                distadd[(size_t)col * m_gr->nbnode + v2] = distadd[(size_t)col * m_gr->nbnode + v] + m_gr->add[(size_t)col * total + i];
              }
              Q.push(std::make_pair(v2, Distances[v2]));
              Visited1[v2]=1;

              touched.push_back(v2);

            }
          }
        }

      }

      if (!Qr.empty()){
        int vv = Qr.top().first;
        double ww = Qr.top().second;
        Qr.pop();

        touched.push_back(vv);

        Visited2[vv]=1;

        if ((Visited1[vv]== 1) && (Distances[vv]+Distances2[vv]) < mu){

          mu=Distances[vv]+Distances2[vv];
          for (int col = 0; col < n_add; col++){
            mu2[col] = distadd[(size_t)col * m_gr->nbnode + vv] + distadd2[(size_t)col * m_gr->nbnode + vv];
          }

        }

        if (ww <= Distances2[vv]) {
          if (Stall_par(vv, Distances2, m_gr->nodeG, m_gr->wG, m_gr->indG)) continue;

          for (int i=m_gr->indGr[vv]; i< m_gr->indGr[vv+1]; i++){
            int vv2 = m_gr->nodeGr[i];
            double ww2 = m_gr->wGr[i];

            if (Distances2[vv] + ww2 < Distances2[vv2]) {
              Distances2[vv2] = Distances2[vv] + ww2;
              for (int col = 0; col < n_add; col++){
                distadd2[(size_t)col * m_gr->nbnode + vv2] = distadd2[(size_t)col * m_gr->nbnode + vv] + m_gr->addr[(size_t)col * totalr + i];
              }
              Qr.push(std::make_pair(vv2, Distances2[vv2]));
              Visited2[vv2]=1;

              touched.push_back(vv2);

            }
          }
        }

      }

    }


    for (int col = 0; col < n_add; col++){
      m_result[(size_t)col * n_pairs + k] = mu2[col];
    }

    for (int i=0; i<(int)touched.size();i++) {
      int node = touched[i];
      Distances[node] = std::numeric_limits<double>::max();
      Distances2[node] = std::numeric_limits<double>::max();
      for (int col = 0; col < n_add; col++){
        distadd[(size_t)col * m_gr->nbnode + node] = std::numeric_limits<double>::max();
        distadd2[(size_t)col * m_gr->nbnode + node] = std::numeric_limits<double>::max();
      }
      Visited1[node] = 0;
      Visited2[node] = 0;
    }

    touched.clear();

    mu=std::numeric_limits<double>::max();
    std::fill(mu2.begin(), mu2.end(), std::numeric_limits<double>::max());
  }

}


