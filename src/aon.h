#include "graph.h"
#include "cgraph.h"
#include <string>
#include <RcppParallel.h>
#include <Rcpp.h>
using namespace std;
using namespace RcppParallel;

struct aonGraph : public RcppParallel::Worker{
  //input
  Graph* m_gr;
  IVec m_dep;
  IVec m_arr;
  DVec m_demand;
  int algorithm;
  //output
  vector<double> m_result;
  vector<pair<int, pair<vector<int>, vector<double>>>> m_od;
  int n_iter;

  // constructor
  aonGraph(Graph* gr, IVec dep, IVec arr, DVec demand, int algo);

  aonGraph(aonGraph &pn, RcppParallel::Split);





  // algorithms
  void dijkstra(std::size_t begin, std::size_t end);
  void dijkstra_reverse(std::size_t begin, std::size_t end);
  void bidir(std::size_t begin, std::size_t end);
  void nba(std::size_t begin, std::size_t end);

  void operator()(std::size_t begin, std::size_t end);


  void join(aonGraph &pb);
};


struct aonGraphC : public RcppParallel::Worker{
  //input
  CGraph* m_gr;
  Graph* m_or;
  IVec m_dep;
  IVec m_arr;
  DVec m_demand;
  int algorithm;
  //output
  vector<double> m_result;
  vector<pair<int, pair<vector<int>, vector<double>>>> m_od;
  int n_iter;

  // constructor
  aonGraphC(CGraph* gr, Graph* original_graph, IVec dep, IVec arr, DVec demand, int algo);

  aonGraphC(aonGraphC &pn, RcppParallel::Split);

  // algorithms
  void dijkstra_bidir(std::size_t begin, std::size_t end);

  void phast(std::size_t begin, std::size_t end);

  // reverse
  void phastr(std::size_t begin, std::size_t end);

  void operator()(std::size_t begin, std::size_t end);


  void join(aonGraphC &pb);
};


// disaggregate aon from contracted to normal graph
struct unpackC : public RcppParallel::Worker{
  CGraph* m_gr; // just for use unpack method and shortcuts
  Graph* m_or; // normal graph
  Graph* m_cg;  // full contracted graph
  DVec m_input; // contracted aon output from aonGraphC
  bool phast_rename;
  DVec m_result; //normal graph aon

  IVec node_dict;

  unpackC(CGraph* gr, Graph* original, Graph* cg, DVec input, bool convert);
  unpackC(unpackC &pn,  RcppParallel::Split);
  // iterate over m_cg->indG
  void operator()(std::size_t begin, std::size_t end);
  void join(unpackC &pb);


};
