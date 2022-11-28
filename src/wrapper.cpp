#include "graph.h"
#include "cgraph.h"
#include "aon.h"
#include "aggc.h"
#include <string>
#include <vector>
#include <RcppParallel.h>
#include <Rcpp.h>
using namespace std;


template <typename Out>
void split(std::string &s, char delim, Out result) {
  std::istringstream iss(s);
  std::string item;
  while (std::getline(iss, item, delim)) {
    *result++ = item;
  }
}

std::vector<std::string> split(std::string &s, char delim) {
  std::vector<std::string> elems;
  split(s, delim, std::back_inserter(elems));
  return elems;
}


void check_nas_vec(Rcpp::NumericVector &vec){
  for (unsigned int i = 0; i < vec.size(); i++){
    if (vec[i] == numeric_limits<double>::max()) vec[i] = Rcpp::NumericVector::get_na();
  }
}

void check_nas_mat(Rcpp::NumericMatrix &mat){
  for (unsigned int i = 0; i < mat.size(); i++){
    if (mat[i] == numeric_limits<double>::max()) mat[i] = Rcpp::NumericVector::get_na();
  }
}

// [[Rcpp::export]]

Rcpp::NumericVector cppdist(std::vector<int> &gfrom, std::vector<int> &gto, std::vector<double> &gw, int nb,
                            std::vector<double> lat, std::vector<double> lon, double k,
                            std::vector<int> dep, std::vector<int> arr, int algo){
  Graph network(gfrom, gto, gw, nb);

  network.to_adj_list(false);

  if (algo == 1 || algo == 3) {
    network.setReverse();
    network.to_adj_list(true);
  }
  if (algo == 2 || algo == 3){
    network.setLatLon(lat, lon);
    network.k = k;
  }

  Rcpp::NumericVector result = Rcpp::wrap(network.routing_dvec(dep, arr, algo));
  check_nas_vec(result);

  return result;
}

// [[Rcpp::export]]
Rcpp::List cpppath(std::vector<int> &gfrom, std::vector<int> &gto, std::vector<double> &gw, int nb,
                            std::vector<double> lat, std::vector<double> lon, double k,
                            std::vector<std::string> dict, std::vector<int> keep,
                            std::vector<int> dep, std::vector<int> arr, double lim, int algo){
  Graph network(gfrom, gto, gw, nb);
  network.to_adj_list(false);
  network.setDict(dict);

  if (algo == 1 || algo == 3 || algo == 5) {
    network.setReverse();
    network.to_adj_list(true);
  }
  if (algo == 2 || algo == 3){
    network.setLatLon(lat, lon);
    network.k = k;
  }

  return Rcpp::wrap(network.routing_svec(dep, arr,keep, lim, algo));
}

// [[Rcpp::export]]
Rcpp::NumericMatrix cppdistmat(std::vector<int> &gfrom, std::vector<int> &gto, std::vector<double> &gw, int nb,
                               std::vector<int> dep, std::vector<int> arr){
  Graph network(gfrom, gto, gw, nb);
  network.to_adj_list(false);

  Rcpp::NumericMatrix result = Rcpp::wrap(network.routing_dmat(dep, arr));
  check_nas_mat(result);
  return result;

}

// [[Rcpp::export]]
Rcpp::List cpppathmat(std::vector<int> &gfrom, std::vector<int> &gto, std::vector<double> &gw, int nb,
                      std::vector<std::string> dict, std::vector<int> keep,
                      std::vector<int> dep, std::vector<int> arr, std::vector<double> lim, bool setdif, int algo, bool reverse){

  Graph network(gfrom, gto, gw, nb);
  network.to_adj_list(false);
  network.setDict(dict);

  vector<SVec> final = network.routing_smat(dep, arr, keep, lim, setdif, algo);

  if (reverse){

    Rcpp::List result(final[0].size());

    for (int i = 0; i < final[0].size(); i++){

      Rcpp::List result2(final.size());

      for (int j = 0; j < final.size(); j++){

        SVec x = split(final[j][i], ' ');
        if (x.size() == 1 && algo == 0){
          SVec s;
          result2[j] = s;
        } else{
          result2[j] = x;
        }


      }
      //SVec().swap(final[i]);
      result[i] = result2;
    }
    return Rcpp::wrap(result);


  } else{
    Rcpp::List result(final.size());
    for (int i = 0; i < final.size(); i++){
      Rcpp::List result2(final[i].size());
      for (int j = 0; j < final[i].size(); j++){

        SVec x = split(final[i][j], ' ');
        std::reverse(x.begin(), x.end());

        if (x.size() == 1 && algo == 0){
          SVec s;
          result2[j] = s;
        } else{
          result2[j] = x;
        }

      }
      SVec().swap(final[i]);
      result[i] = result2;
    }
    return Rcpp::wrap(result);
  }


}


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::List cppsimplify(std::vector<int> &gfrom, std::vector<int> &gto, std::vector<double> &gw, int nb,
                       std::vector<int> keep, bool rm_loop, bool iterate, bool progress){

  Graph network(gfrom, gto, gw, nb, true);

  network.simplify(rm_loop, keep, iterate, progress);
  Rcpp::List res = network.getEdges();

  return res;
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::List cppcontract(std::vector<int> &gfrom,std::vector<int> &gto,std::vector<double> &gw,int NbNodes,bool display_progress){
  CGraph network(gfrom, gto, gw, NbNodes);

  network.contract(display_progress);
  Rcpp::List final(3);
  Rcpp::List res = network.getEdges();
  Rcpp::List sc(3);
  sc[0] = network.shortf;
  sc[1] = network.shortt;
  sc[2] = network.shortc;
  final[0] = res;
  final[1] = network.rank;
  final[2] = sc;

  return final;
}


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::NumericVector cppdistC(std::vector<int> &gfrom, std::vector<int> &gto, std::vector<double> &gw, int nb,
                             std::vector<int> &rank, std::vector<int> &shortf, std::vector<int> &shortt, std::vector<int> &shortc,
                            bool phast,
                            std::vector<int> dep, std::vector<int> arr, int algo){
  CGraph network(gfrom, gto, gw, nb, rank, shortf, shortt, shortc, phast);
  network.to_adj_list(false, false);
  network.to_adj_list(true, false);

  Rcpp::NumericVector result = Rcpp::wrap(network.routing_dvec(dep, arr, algo));
  check_nas_vec(result);

  return result;
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::List cpppathC(std::vector<int> &gfrom, std::vector<int> &gto, std::vector<double> &gw, int nb,
                             std::vector<int> &rank, std::vector<int> &shortf, std::vector<int> &shortt, std::vector<int> &shortc,
                             bool phast,
                             std::vector<std::string> dict, std::vector<int> keep,
                             std::vector<int> dep, std::vector<int> arr, int algo){
  CGraph network(gfrom, gto, gw, nb, rank, shortf, shortt, shortc, phast);
  network.setDict(dict);
  network.construct_shortcuts();
  network.to_adj_list(false, false);
  network.to_adj_list(true, false);

  return Rcpp::wrap(network.routing_svec(dep, arr, keep, algo));
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::NumericMatrix cppdistmatC(std::vector<int> &gfrom, std::vector<int> &gto, std::vector<double> &gw, int nb,
                                std::vector<int> &rank, std::vector<int> &shortf, std::vector<int> &shortt, std::vector<int> &shortc,
                                bool phast,
                                std::vector<int> dep, std::vector<int> arr, int algo){

  CGraph network(gfrom, gto, gw, nb, rank, shortf, shortt, shortc, phast);

  Rcpp::NumericMatrix result = Rcpp::wrap(network.routing_dmat(dep, arr, algo));
  check_nas_mat(result);
  return result;
}


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::List cpppathmatC(std::vector<int> &gfrom, std::vector<int> &gto, std::vector<double> &gw, int nb,
                                std::vector<int> &rank, std::vector<int> &shortf, std::vector<int> &shortt, std::vector<int> &shortc,
                                bool phast,
                                std::vector<std::string> dict, std::vector<int> keep,
                                std::vector<int> dep, std::vector<int> arr){

  CGraph network(gfrom, gto, gw, nb, rank, shortf, shortt, shortc, phast);
  network.setDict(dict);
  network.construct_shortcuts();
  network.to_adj_list(false, true);
  network.to_adj_list(true, true);

  vector<SVec> final = network.routing_smat(dep, arr, keep);

  Rcpp::List result(final.size());

  for (int i = 0; i < final.size(); i++){
    Rcpp::List result2(final[i].size());
    for (int j = 0; j < final[i].size(); j++){

      SVec x = split(final[i][j], ' ');
      result2[j] = x;

    }
    SVec().swap(final[i]);
    result[i] = result2;
  }

  return result;
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::NumericMatrix cpppadd(std::vector<int> &gfrom, std::vector<int> &gto, std::vector<double> &gw, std::vector<double> &gadd, int nb,
            std::vector<int> dep, std::vector<int> arr){
  Graph network(gfrom, gto, gw, gadd, nb);
  Rcpp::NumericMatrix result = Rcpp::wrap(network.routing_dmat(dep, arr));
  check_nas_mat(result);
  return result;
}


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::NumericVector cppdistadd(std::vector<int> &gfrom, std::vector<int> &gto, std::vector<double> &gw, std::vector<double> &gadd,int nb,
                            std::vector<double> lat, std::vector<double> lon, double k,
                            std::vector<int> dep, std::vector<int> arr, int algo){
  Graph network(gfrom, gto, gw, gadd, nb);


  if (algo == 2 || algo == 3){
    network.setLatLon(lat, lon);
    network.k = k;
  }

  Rcpp::NumericVector result = Rcpp::wrap(network.routing_dvec(dep, arr, algo));
  check_nas_vec(result);

  return result;
}


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::NumericVector cppdistaddC(std::vector<int> &orfrom, std::vector<int> &orto, std::vector<double> &orw, std::vector<double> &gadd,
                                std::vector<int> &gfrom, std::vector<int> &gto, std::vector<double> &gw, int nb,
                                std::vector<int> &rank, std::vector<int> &shortf, std::vector<int> &shortt, std::vector<int> &shortc,
                                bool phast,
                                std::vector<int> dep, std::vector<int> arr, int algo){
  Graph network(orfrom, orto, orw, gadd, nb);

  CGraph networkc(gfrom, gto, gw, nb, rank, shortf, shortt, shortc, phast);
  networkc.construct_shortcuts();
  networkc.to_adj_list(false, phast);
  networkc.to_adj_list(true, phast);

  aggC dijfunc(&networkc, &network);
  parallelFor(0, networkc.nbnode, dijfunc);

  networkc.add = dijfunc.m_result;
  networkc.addr = dijfunc.m_result2;

  Rcpp::NumericVector result = Rcpp::wrap(networkc.routing_dvec(dep, arr, algo));
  check_nas_vec(result);

  return result;
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::NumericMatrix cppaddC(std::vector<int> &orfrom, std::vector<int> &orto, std::vector<double> &orw, std::vector<double> &gadd,
                       std::vector<int> &gfrom, std::vector<int> &gto, std::vector<double> &gw, int nb,
                       std::vector<int> &rank, std::vector<int> &shortf, std::vector<int> &shortt, std::vector<int> &shortc,
                       bool phast,
                       std::vector<int> dep, std::vector<int> arr, int algo){
  Graph network(orfrom, orto, orw, gadd, nb);

  CGraph networkc(gfrom, gto, gw, nb, rank, shortf, shortt, shortc, phast);
  networkc.construct_shortcuts();
  networkc.to_adj_list(false, phast);
  networkc.to_adj_list(true, phast);

  aggC dijfunc(&networkc, &network);
  parallelFor(0, networkc.nbnode, dijfunc);

  networkc.add = dijfunc.m_result;
  networkc.addr = dijfunc.m_result2;

  Rcpp::NumericMatrix result = Rcpp::wrap(networkc.routing_dmat(dep, arr, algo));
  check_nas_mat(result);

  return result;
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::List cppaon(std::vector<int> &gfrom, std::vector<int> &gto, std::vector<double> &gw,int nb,
                           std::vector<double> lat, std::vector<double> lon, double k,
                           std::vector<int> dep, std::vector<int> arr, std::vector<double> dem,
                           int algo){
  Graph network(gfrom, gto, gw, nb);
  network.to_adj_list(false);
  if (algo > 0) {
    network.setReverse();
    network.to_adj_list(true);
  }
  if (algo == 3){
    network.setLatLon(lat, lon);
    network.k = k;
  }

  Rcpp::NumericVector result = Rcpp::wrap(network.getaon(dep, arr, dem, k, algo));
  check_nas_vec(result);

  Rcpp::List final(4);

  IVec from(network.nbedge);
  int count = 0;
  for (int i = 0; i < (network.indG.size() - 1); i++){
    for (int j = network.indG[i]; j < network.indG[i+1]; j++){
      from[count] = i;
      count++;
    }
  }

  final[0] = Rcpp::wrap(from);
  final[1] = Rcpp::wrap(network.nodeG);
  final[2] = Rcpp::wrap(network.wG);
  final[3] = result;

  return final;
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::List cppaonC(std::vector<int> &orfrom, std::vector<int> &orto, std::vector<double> &orw,
                            std::vector<int> &gfrom, std::vector<int> &gto, std::vector<double> &gw, int nb,
                            std::vector<int> &rank, std::vector<int> &shortf, std::vector<int> &shortt, std::vector<int> &shortc,
                            bool phast,
                           std::vector<int> dep, std::vector<int> arr, std::vector<double> dem,
                           int algo){

  Graph original(orfrom, orto, orw, nb);
  original.to_adj_list(false);

  Graph contracted(gfrom, gto, gw, nb);
  contracted.to_adj_list(false);

  CGraph network(gfrom, gto, gw, nb, rank, shortf, shortt, shortc, phast);
  network.construct_shortcuts();
  network.to_adj_list(false, phast);
  network.to_adj_list(true, phast);

  Rcpp::List result(4);
  IVec from(original.nbedge);
  int count = 0;
  for (int i = 0; i < (original.indG.size() - 1); i++){
    for (int j = original.indG[i]; j < original.indG[i+1]; j++){
      from[count] = i;
      count++;
    }
  }

  DVec flow = network.getaon(&contracted, dep, arr, dem, algo);


  unpackC dijfunc(&network, &original, &contracted, flow, false); // false, because converted from R side
  parallelReduce(0, contracted.indG.size()-1, dijfunc);

  result[0] = Rcpp::wrap(from);
  result[1] = Rcpp::wrap(original.nodeG);
  result[2] = Rcpp::wrap(original.wG);
  result[3] = Rcpp::wrap(dijfunc.m_result);

  return result;
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::List cpptraffic(std::vector<int> &gfrom, std::vector<int> &gto, std::vector<double> &gw,
                      std::vector<double> &gflow,  std::vector<double> &gaux,  std::vector<double> &gftt,
                      std::vector<double> &galpha,  std::vector<double> &gbeta,  std::vector<double> &gcap,
                      int nb,
                      std::vector<double> lat, std::vector<double> lon, double k,
                      std::vector<int> dep, std::vector<int> arr, std::vector<double> dem,
                      double max_gap, int max_it, int method, int aon_method,
                      bool contract, bool phast, bool verbose){

  Graph network(gfrom,gto, gw, gflow,gaux, gftt, galpha, gbeta, gcap, nb);
  network.setLatLon(lat, lon);
  network.k= k;
  network.assign_traffic(method, max_gap, max_it, dep, arr, dem, aon_method, contract, phast, verbose);

  Rcpp::List result(10);
  result[0] = network.indG2;
  result[1] = network.nodeG;
  result[2] = network.ftt;
  result[3] = network.wG;
  result[4] = network.flow;
  result[5] = network.cap;
  result[6] = network.alpha;
  result[7] = network.beta;
  result[8] = network.gap;
  result[9] = network.it;
  return (result);

}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::List cppalgB(std::vector<int> &gfrom, std::vector<int> &gto, std::vector<double> &gw,
                      std::vector<double> &gflow,  std::vector<double> &gaux,  std::vector<double> &gftt,
                      std::vector<double> &galpha,  std::vector<double> &gbeta,  std::vector<double> &gcap,
                      int nb,
                      std::vector<double> lat, std::vector<double> lon, double k,
                      std::vector<int> dep, std::vector<int> arr, std::vector<double> dem,
                      double max_gap, int max_it, int aon_method,
                      int batch_size, int n_batch, std::string file_path,
                      int inner_iter,
                      double NUM_TOL,
                      bool contract,
                      bool phast,
                      bool verbose){

  Graph network(gfrom,gto, gw, gflow,gaux, gftt, galpha, gbeta, gcap, nb);
  network.setLatLon(lat, lon);
  network.k= k;
  network.setReverse();
  network.to_adj_list(true);

  network.algorithmB(batch_size, n_batch, file_path, max_gap, max_it, dep, arr, dem, aon_method, inner_iter, NUM_TOL, contract, phast, verbose );

  Rcpp::List result(10);
  result[0] = network.indG2;
  result[1] = network.nodeG;
  result[2] = network.ftt;
  result[3] = network.wG;
  result[4] = network.flow;
  result[5] = network.cap;
  result[6] = network.alpha;
  result[7] = network.beta;
  result[8] = network.gap;
  result[9] = network.it;

  return (result);

}
