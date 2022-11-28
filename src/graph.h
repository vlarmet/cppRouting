#ifndef GRAPH_H
#define GRAPH_H

#include <queue>
#include <string>
#include <RcppParallel.h>
#include <Rcpp.h>
using namespace std;

// Priority queue comparator
struct comp{

  bool operator()(const pair<int, double> &a, const pair<int, double> &b){
    return a.second > b.second;
  }
};

// Data structures
using G = vector<vector<pair<int, double> > >;
using PQ = priority_queue<pair<int, double>, vector<pair<int, double> >, comp >;
using IVec = vector<int>;
using DVec = vector<double>;
using SVec = vector<string>;
using VVec = vector<SVec >;
using VIVec = vector<IVec >;

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]



// Graph object

class Graph{
public:
  // Mandatory attribute
  int nbnode;
  int nbedge;
  G data;
  int n_core;

  // Optionnal
  DVec lat;
  DVec lon;
  double k;
  SVec dict;
  G dataR;

  // Parallel graph : 3 vectors (adjacency list)
  IVec nodeG;
  IVec indG;
  DVec wG;
  IVec indG2; // used to update graphs weights (traffic assignment)
  // edge description for traffic assignment
  DVec flow;
  DVec aux;
  DVec ftt;
  DVec alpha;
  DVec beta;
  DVec cap;

  // additional weight
  DVec add;
  DVec addr;

  /// reversed
  IVec nodeGr;
  IVec indGr;
  DVec wGr;

  int it;
  double gap;
  double tstt;
  double sptt;

  // Constructor
  Graph(IVec &gfrom, IVec &gto, DVec &gw, int nb);

  // Constructor + bool for removing duplicated edges (useful for simplify)
  Graph(IVec &gfrom, IVec &gto, DVec &gw, int nb, bool rm_dup);

  // Constructor with additional weight
  Graph(IVec &gfrom, IVec &gto, DVec &gw, DVec &gadd, int nb);

  // constructor with full edges description
  Graph(IVec &gfrom, IVec &gto, DVec &gw,
        DVec &gflow, DVec &gaux, DVec &gftt,
        DVec &galpha, DVec &gbeta, DVec &gcap, int nb);

  // Setters
  void setReverse();
  void setLatLon(DVec &m_lat, DVec &m_lon);
  void setDict(SVec &m_dict);
  void to_adj_list(bool reversed);

  // return edge description from current graph
  Rcpp::List getEdges();

  // routing algorithms, grouped by return type
  DVec routing_dvec(IVec dep, IVec arr, int algo);

  Rcpp::NumericMatrix routing_dmat(IVec dep, IVec arr);

  // for multi path and multi isochrone, strings are concatened to avoid huge memory fragmentation
  vector<SVec> routing_smat(IVec dep, IVec arr, IVec keep, DVec lim, bool setdif, int algo);

  vector<SVec> routing_svec(IVec dep, IVec arr, IVec keep, double lim, int algo);


  // Simplification
  void simp(vector<IVec> &Gr,
            IVec &Junction,
            vector<IVec > &Edges,
            IVec &keep,
            IVec &Treated,
            bool &loop);

  void simplify(bool loop,IVec keep,bool iterate, bool progress);

  // Traffic assignment
  DVec getaon(IVec dep, IVec arr, DVec dem, double k_, int mode);


  double f(double theta);
  double bissection(double tol);
  void update_flow(double theta);
  void update_cost();
 // void update_cost2(); // update cost in data attributes from adjacency list (useful when CH is used in traffic assignment)
  void update_tstt();
  void update_sptt();
  void cfw_update_aux(DVec &aux_m1, DVec &d);
  void bfw_update_aux(DVec &aux_m1, DVec &aux_m2, double theta_old, DVec &d, DVec &d2);

  void assign_traffic(int method, double max_gap, int max_iter,
                      IVec dep, IVec arr, DVec dem, int mode,
                      bool contract, bool phast, bool verbose);

  void algorithmB(int batch_size, int n_batch, string file_path, double max_gap, int max_iter,
                         IVec dep, IVec arr, DVec dem, int mode, int inner_iter, double NUM_TOL,
                         bool contract, bool phast, bool verbose);

};



#endif
