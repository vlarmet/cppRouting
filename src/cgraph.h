#ifndef CGRAPH_H
#define CGRAPH_H

#include "graph.h"
#include <vector>
#include <RcppParallel.h>
#include <Rcpp.h>
using namespace std;


class CGraph{
public:

  int nbnode;
  int nbedge;
  G cdata; // contracted graph
  G cdataR;
  bool is_contracted;
  IVec rank;
  IVec shortf;
  IVec shortt;
  IVec shortc;
  vector<vector<pair<int, int >>> shortcuts;
  SVec dict;

  // Parallel graph : 3 vectors (adjacency list)
  IVec nodeG;
  IVec indG;
  DVec wG;

  /// reversed
  IVec nodeGr;
  IVec indGr;
  DVec wGr;

  // additional weight
  DVec add;
  DVec addr;

  // constructors
  CGraph(Graph* graph); //from a c++ Graph object
  // from R
  CGraph(IVec &gfrom, IVec &gto, DVec &gw, int nb, IVec &m_rank, IVec &m_shortf, IVec &m_shortt, IVec &m_shortc, bool phast); // from R (already contracted)
  CGraph(IVec &gfrom, IVec &gto, DVec &gw, int nb);

  G getReverse(G &graph);
  void setDict(SVec &m_dict);
  Rcpp::List getEdges();
  void construct_shortcuts();
  void unpack(IVec &path);
  void to_adj_list(bool reversed, bool sorting);


  // contraction hierarchies functions
  /// remove an edge if there is a shorter path between "from" and "to" nodes (used for initially cleaning the graph)
  void clean(G &Graph,
             G &OrGraph,
             int dep, IVec &arr, DVec &lim,
             DVec &distances);
  /// return the number of shortcuts created if a given node is contracted
  int find_shortcuts(G &graph,
                     int dep, IVec &arr,
                     DVec &lim, int node,
                     DVec &distances);
  /// Heuristic used for ordering :  nb of shortcuts potentially created - nb of incoming edge - nb of outcoming edge
  int edge_dif(int node, G &graph, G &graphr,
               DVec &distances);
  /// Store created shortcuts
  vector<pair<int, pair<int,double> > > get_shortcuts(G &graph,
                                                      int dep,
                                                      IVec &arr,
                                                      DVec &lim,
                                                      int node,
                                                      DVec &distances,
                                                      bool reversed);
  /// Contract one node : create shortcuts, remove node from current graph, add shortcuts ...
  void contract_one_node(int node,
                         G &graph,
                         G &graphr,
                         G &cgraph,
                         DVec &distances,
                         IVec &contracted);

  /// main function, contract the whole network
  void contract(bool display_progress);

  // Stall on demand (reduce search space during modified dijkstra algorithm)
  bool stall(int &node, DVec &distances, G &graph);


  // routing algorithms
  DVec routing_dvec(IVec dep, IVec arr, int algo);
  vector<SVec> routing_svec(IVec dep, IVec arr, IVec keep, int algo);
  Rcpp::NumericMatrix routing_dmat(IVec dep, IVec arr, int algo);
  vector<SVec> routing_smat(IVec dep, IVec arr, IVec keep);

  // All or Nothing assignment
  DVec getaon(Graph* original_graph, IVec dep, IVec arr, DVec dem, int algo);

};



#endif
