#include "graph.h"
#include "cgraph.h"
#include "aon.h"
#include "bush.h"
#include <fstream>
#include <string>
#include <chrono>
#include <RcppParallel.h>
#include <Rcpp.h>
using namespace std;
using namespace RcppParallel;


double Graph::bissection(double tol){
  double inf = 0.0;
  double sup = 1.0;
  double tentative = -1.0;
  if ((f(inf) * f(sup)) < 0){
    while ((sup - inf) > tol ){
      tentative = (inf + sup)/2;
      if (f(tentative) > 0){
        sup = tentative;
      } else{
        inf = tentative;
      }
    }
  }


  return tentative;

}

double Graph::f(double theta){
  double result = 0;

  for (int i = 0; i < flow.size(); i++){
    result += (ftt[i] * (1 + alpha[i]*(pow(((theta * aux[i]) + ((1 - theta) * flow[i]))/cap[i], beta[i])))) * (aux[i] - flow[i]);
  }

  return result;
}

void Graph::update_flow(double theta){
  for (int i = 0; i < nbedge; i++){
    double new_flow = theta * aux[i] + (1.0 - theta) * flow[i];
    flow[i] = new_flow;

  }
}

void Graph::update_cost(){
  for (int i = 0; i < nbedge; i++){
    // adjacency list
    double old_cost = wG[i];
    double new_cost = ftt[i] * (1.0 + alpha[i] * (pow(flow[i]/cap[i], beta[i])));
    wG[i] = new_cost;
    // reversed
    for (int j = indGr[nodeG[i]]; j < indGr[nodeG[i] + 1]; j++){
      if (nodeGr[j] == indG2[i] && wGr[j] == old_cost){
        wGr[j] = new_cost;
        break;
      }
    }

    // data structure
    for (int j = 0; j < data[indG2[i]].size(); j++){
      if (data[indG2[i]][j].first == nodeG[i] && data[indG2[i]][j].second == old_cost){
        data[indG2[i]][j].second = new_cost;
        break;
      }
    }

  }
}



void Graph::update_tstt(){
  double result = 0;
  for (int i = 0; i < nbedge; i++){
    result += flow[i] * wG[i];
  }

  tstt = result;
}

void Graph::update_sptt(){
  double result = 0;
  for (int i = 0; i < nbedge; i++){
    result += aux[i] * wG[i];
  }

  sptt = result;
}

void Graph::cfw_update_aux(DVec &aux_m1, DVec &d){
  // alpha determination

  // DVec d(nbedge);
  for (int i = 0; i < nbedge; i++){
    d[i] = ftt[i] * (alpha[i] * beta[i] * pow(flow[i]/cap[i], beta[i] - 1) )/cap[i];
  }

  double num = 0.0;
  double den = 0.0;

  for (int i = 0; i < nbedge; i++){
    num += (aux_m1[i] - flow[i]) * d[i] * (aux[i] - flow[i]);
    den += (aux_m1[i] - flow[i]) * d[i] * ((aux[i] - flow[i]) - (aux_m1[i] - flow[i]));
  }
  //  double test = ftt[0] * (alpha[0] * beta[0] * pow(flow[0]/cap[0], beta[0] - 1) );
  double alpha = num/(den + numeric_limits<double>::epsilon());
  //  Rcpp::Rcout << "alpha :" << alpha << " den:" << den << " d:" << d[0] << " cap:"<< cap[0] << " flow:"<< flow[0] << " numd:"<< test <<" pow:"<<pow(flow[0]/cap[0], beta[0] - 1)<<endl;
  if (alpha > 0.99) alpha = 0.99;
  if (alpha < 0) alpha = 0;


  // AON update

  for (int i = 0; i < nbedge; i++){
    aux[i] = alpha * aux_m1[i] + (1.0 - alpha) * aux[i];
  }
}

void Graph::bfw_update_aux(DVec &aux_m1, DVec &aux_m2, double theta_old, DVec &d, DVec &d2){

  //  DVec d(nbedge);
  //  DVec d2(nbedge);

  for (int i = 0; i < nbedge; i++){
    d[i] = ftt[i] * (alpha[i] * beta[i] * pow(flow[i]/cap[i], beta[i] - 1) )/cap[i];
    if (flow[i] == 0) d[i] = ftt[i];
    d2[i] = theta_old * aux_m1[1] - flow[i] + (1.0 - theta_old) * aux_m2[i];
  }

  // mu
  double mu_num = 0;
  double mu_den = 0;

  for (int i = 0; i < nbedge; i++){
    mu_num += d2[i] * d[i] * (aux[i] - flow[i]);
    mu_den += d2[i] * d[i] * (aux_m2[i] - aux_m1[i]);
  }
  double mu = -mu_num/mu_den;

  // nu
  double nu_num = 0;
  double nu_den = 0;

  for (int i = 0; i < nbedge; i++){
    nu_num += (aux_m1[i] - flow[i]) * d[i] * (aux[i] - flow[i]);
    nu_den += (aux_m1[i] - flow[i]) * d[i] * (aux_m1[i] - flow[i]);
  }
  double nu = -nu_num/nu_den + ((mu * theta_old)/(1.0 - theta_old));

  // beta coefficients
  double beta0 = 1/(1.0 + mu + nu);
  double beta1 = beta0 * nu;
  double beta2 = beta0 * mu;

  if (beta0 < 0 || beta1 < 0 || beta2 < 0){
    beta0 = 1.0;
    beta1 = 0.0;
    beta2 = 0.0;
  }

  // AON update
  for (int i = 0; i < nbedge; i++){
    aux[i] = beta0 * aux[i] + beta1 * aux_m1[i] + beta2 * aux_m2[i];
  }

}

DVec Graph::getaon(IVec dep, IVec arr, DVec dem, double k_, int mode){
  k = k_;
  aonGraph dijfunc(this, dep, arr, dem, mode);
  parallelReduce(0, dijfunc.n_iter, dijfunc);
  return dijfunc.m_result;
}


DVec CGraph::getaon(Graph* original_graph, IVec dep, IVec arr, DVec dem, int algo){
  aonGraphC dijfunc(this, original_graph, dep, arr, dem, algo);
  parallelReduce(0, dijfunc.n_iter, dijfunc);
  return dijfunc.m_result;
}

// contract graph and compute AON from scratch
DVec ch_aon(Graph* ptr, IVec dep, IVec arr, DVec &dem, bool phast, int algo){
  CGraph cgraph(ptr);
  cgraph.contract(false);
  Rcpp::List cedges = cgraph.getEdges();
  IVec gfrom = cedges[0];
  IVec gto = cedges[1];
  DVec gw = cedges[2];
  // rename node according to its rank
  if (phast){
    for (int i = 0; i < gfrom.size(); i++){
      gfrom[i] = cgraph.nbnode - cgraph.rank[gfrom[i]];
      gto[i] = cgraph.nbnode - cgraph.rank[gto[i]];
    }

    for (int i = 0; i < cgraph.shortf.size(); i++){
      cgraph.shortf[i] = cgraph.nbnode - cgraph.rank[cgraph.shortf[i]];
      cgraph.shortt[i] = cgraph.nbnode - cgraph.rank[cgraph.shortt[i]];
      cgraph.shortc[i] = cgraph.nbnode - cgraph.rank[cgraph.shortc[i]];
    }

    for (int i = 0; i < dep.size(); i++){
      dep[i] = cgraph.nbnode - cgraph.rank[dep[i]];
      arr[i] = cgraph.nbnode - cgraph.rank[arr[i]];
    }
  }
  Graph aug_graph(gfrom, gto, gw, cgraph.nbnode);
  aug_graph.to_adj_list(false);
  CGraph routing_cgraph(gfrom, gto, gw, cgraph.nbnode, cgraph.rank, cgraph.shortf, cgraph.shortt, cgraph.shortc, phast);
  routing_cgraph.construct_shortcuts();
  routing_cgraph.to_adj_list(false, phast);
  routing_cgraph.to_adj_list(true, phast);
  DVec contracted_aux = routing_cgraph.getaon(&aug_graph, dep, arr, dem, algo);

  unpackC dijfunc(&routing_cgraph, ptr, &aug_graph, contracted_aux, phast);
  parallelReduce(0, aug_graph.indG.size()-1, dijfunc);
  return dijfunc.m_result;
}


void Graph::assign_traffic(int method, double max_gap, int max_iter,
                           IVec dep, IVec arr, DVec dem, int mode,
                           bool contract, bool phast,
                           bool verbose){ // 0 : msa, 1 : fw, 2 : cfw, 3 : bfw
  setReverse();
  to_adj_list(true);

  aonGraph gr(this, dep, arr, dem, mode);

  // MSA
  if (method == 0)
  { gap = numeric_limits<double>::max();

    if (contract) {
      aux = ch_aon(this, dep, arr, dem, phast, mode);
    } else{
      parallelReduce(0, gr.n_iter, gr);
      aux = gr.m_result;
      fill(gr.m_result.begin(), gr.m_result.end(), 0.0);
    }
    it = 1;

    while (gap > max_gap){
      double theta = 1.0/it;

      Rcpp::checkUserInterrupt();

      update_flow(theta);
      update_cost();

      if (contract) {
        aux = ch_aon(this, dep, arr, dem, phast, mode);
      } else{
        parallelReduce(0, gr.n_iter, gr);
        aux = gr.m_result;
        fill(gr.m_result.begin(), gr.m_result.end(), 0.0);
      }

      update_tstt();
      update_sptt();
      gap = abs(tstt/sptt - 1);
      if (verbose) Rcpp::Rcout << "iteration " << it << " : " << gap << endl;
      it += 1;
      if (it > max_iter) break;

    }

  }

  // Frank Wolfe
  if (method == 1)
  { gap = numeric_limits<double>::max();


    if (contract) {
      aux = ch_aon(this, dep, arr, dem, phast, mode);
    } else{
      parallelReduce(0, gr.n_iter, gr);
      aux = gr.m_result;
      fill(gr.m_result.begin(), gr.m_result.end(), 0.0);
    }

    it = 1;
    double theta = 0.0;

    while (gap > max_gap){
      theta = bissection(numeric_limits<double>::epsilon());
      if (theta < 0) theta = 1.0/it;

      Rcpp::checkUserInterrupt();

      update_flow(theta);
      update_cost();

      if (contract) {
        aux = ch_aon(this, dep, arr, dem, phast, mode);
      } else{
        parallelReduce(0, gr.n_iter, gr);
        aux = gr.m_result;
        fill(gr.m_result.begin(), gr.m_result.end(), 0.0);
      }

      update_tstt();
      update_sptt();
      gap = abs(tstt/sptt - 1);
      if (verbose) Rcpp::Rcout << "iteration " << it << " : " << gap << endl;
      it += 1;
      if (it > max_iter) break;

    }
  }

  // Conjugate Frank Wolfe
  if (method == 2)
  {
    gap = numeric_limits<double>::max();

    if (contract) {
      aux = ch_aon(this, dep, arr, dem, phast, mode);
    } else{
      parallelReduce(0, gr.n_iter, gr);
      aux = gr.m_result;
      fill(gr.m_result.begin(), gr.m_result.end(), 0.0);
    }

    it = 1;
    double theta = 0.0;
    DVec aux_prev = aux;
    DVec d(nbedge);

    while (gap > max_gap){
      if (it > 1) cfw_update_aux(aux_prev, d);

      Rcpp::checkUserInterrupt();

      theta = bissection(numeric_limits<double>::epsilon());
      //Rcpp::Rcout << theta << endl;
      if (theta < 0) theta = 1.0/it;

      update_flow(theta);
      update_cost();
      aux_prev = aux;

      if (contract) {
        aux = ch_aon(this, dep, arr, dem, phast, mode);
      } else{
        parallelReduce(0, gr.n_iter, gr);
        aux = gr.m_result;
        fill(gr.m_result.begin(), gr.m_result.end(), 0.0);
      }

      update_tstt();
      update_sptt();

      gap = abs(tstt/sptt - 1);
      if (verbose) Rcpp::Rcout << "iteration " << it << " : " << gap << endl;
      it += 1;
      if (it > max_iter) break;

    }

  }
  // Biconjugate Frank Wolfe
  if (method == 3)
  {
    gap = numeric_limits<double>::max();

    if (contract) {
      aux = ch_aon(this, dep, arr, dem, phast, mode);
    } else{
      parallelReduce(0, gr.n_iter, gr);
      aux = gr.m_result;
      fill(gr.m_result.begin(), gr.m_result.end(), 0.0);
    }


    it = 1;
    double theta = 0.0;
    DVec aux_prev = aux;
    DVec aux_prev2 = aux;
    DVec d(nbedge);
    DVec d2(nbedge);

    while (gap > max_gap){
      if (it < 3) cfw_update_aux(aux_prev, d);
      if (it >= 3) bfw_update_aux(aux_prev, aux_prev2, theta, d, d2);

      Rcpp::checkUserInterrupt();

      theta = bissection(numeric_limits<double>::epsilon());
      if (theta < 0) theta = 1.0/it;

      update_flow(theta);
      update_cost();
      aux_prev2 = aux_prev;
      aux_prev = aux;

      if (contract) {
        aux = ch_aon(this, dep, arr, dem, phast, mode);
      } else{
        parallelReduce(0, gr.n_iter, gr);
        aux = gr.m_result;
        fill(gr.m_result.begin(), gr.m_result.end(), 0.0);
      }

      update_tstt();
      update_sptt();
      gap = abs(tstt/sptt - 1);
      if (verbose) Rcpp::Rcout << "iteration " << it << " : " << gap << endl;
      it += 1;
      if (it > max_iter) break;

    }

  }


}


void storeBushes(vector<Bush> &bush_vec, vector<IVec> &batch_vec, int batch_num, string &file_path, IVec &equ, bool equilibrate){
  std::ofstream ofs;
  const char* pointer;
  size_t bytes;

  for (int i = 0; i < batch_vec[batch_num].size(); i++){
    int bush_index = batch_vec[batch_num][i];

    if (equilibrate && equ[bush_index] == 1) continue;

    // edges
    ofs.open(file_path + "edges" + std::to_string(bush_index), std::ios::out | std::ofstream::binary);

    pointer = reinterpret_cast<const char*>(&bush_vec[i].edges[0]);
    bytes = bush_vec[i].edges.size() * sizeof(bush_vec[i].edges[0]);
    ofs.write(pointer, bytes);
    ofs.close();

    // flow
    ofs.open(file_path + "flow" + std::to_string(bush_index), std::ios::out | std::ofstream::binary);

    pointer = reinterpret_cast<const char*>(&bush_vec[i].flow[0]);
    bytes = bush_vec[i].flow.size() * sizeof(bush_vec[i].flow[0]);
    ofs.write(pointer, bytes);
    ofs.close();

    // order
    ofs.open(file_path + "order" + std::to_string(bush_index), std::ios::out | std::ofstream::binary);

    pointer = reinterpret_cast<const char*>(&bush_vec[i].order[0]);
    bytes = bush_vec[i].order.size() * sizeof(bush_vec[i].order[0]);
    ofs.write(pointer, bytes);
    ofs.close();
  }

}

void loadBushes(vector<Bush> &bush_vec, vector<IVec> &batch_vec, int batch_num, string &file_path, IVec &equ, bool equilibrate){
  std::ifstream ifs;
  char* pointer;
  size_t bytes;

  for (int i = 0; i < batch_vec[batch_num].size(); i++){
    int bush_index = batch_vec[batch_num][i];

    if (equilibrate && equ[bush_index] == 1) continue;

    // edges
    ifs.open(file_path + "edges" + std::to_string(bush_index), std::ios::in | std::ifstream::binary);

    pointer = reinterpret_cast<char*>(&bush_vec[i].edges[0]);
    bytes = bush_vec[i].edges.size() * sizeof(bush_vec[i].edges[0]);
    ifs.read(pointer, bytes);
    ifs.close();

    // flow
    ifs.open(file_path + "flow" + std::to_string(bush_index), std::ios::in | std::ifstream::binary);

    pointer = reinterpret_cast<char*>(&bush_vec[i].flow[0]);
    bytes = bush_vec[i].flow.size() * sizeof(bush_vec[i].flow[0]);
    ifs.read(pointer, bytes);
    ifs.close();

    // order
    ifs.open(file_path + "order" + std::to_string(bush_index), std::ios::in | std::ifstream::binary);

    pointer = reinterpret_cast<char*>(&bush_vec[i].order[0]);
    bytes = bush_vec[i].order.size() * sizeof(bush_vec[i].order[0]);
    ifs.read(pointer, bytes);
    ifs.close();

  }


}


void Graph::algorithmB(int batch_size,int n_batch, string file_path, double max_gap, int max_iter,
                       IVec dep, IVec arr, DVec dem, int mode,
                       int inner_iter,
                       double NUM_TOL,
                       bool contract,
                       bool phast,
                       bool verbose){

  // map representing OD matrix: first : origin, second : pair of (vector) destination nodes and demands
  map<int, pair<IVec, DVec>> od_map;
  map<int, pair<IVec, DVec>>::iterator ite;

  for (int i = 0; i < dep.size(); i++) {
    od_map[dep[i]].first.push_back(arr[i]);
    od_map[dep[i]].second.push_back(dem[i]);

  }

  // convert to a vector
  vector<pair<int, pair<vector<int>, vector<double>>>> od(od_map.size());
  int ind = 0;
  for (ite = od_map.begin(); ite != od_map.end(); ite++){
    od[ind] = make_pair(ite->first, ite->second);
    ind++;
  }

  // Batch vector : store map indexes for each batch
  vector<IVec> batches(n_batch);

  int num_batch = 0;
  int nb_bushes = 0;

  for (int i = 0; i < od.size(); i++){
    batches[num_batch].push_back(i);
    nb_bushes++;
    if (nb_bushes == batch_size){
      num_batch++;
      nb_bushes = 0;
    }
  }


  IVec equilibrated(od.size(), 0);
  // Bushes initialization
  /// shared vectors
  Bush_vectors bv(nbnode);

  /// vector of bushes
  vector<Bush> bushes(batch_size);

  if (verbose) Rcpp::Rcout << "Bushes initialization..." <<  endl;
  for (int i = 0; i < batches.size(); i++){
    for (int j = 0; j < batches[i].size(); j++){

      int index = batches[i][j];

      Bush b(this, od[index].first, od[index].second.first, od[index].second.second, &bv, NUM_TOL);
      b.ordering();
      b.minmaxtree2(1);
      b.loadAON();

      bushes[j] = b;
    }
    if (n_batch > 1) storeBushes(bushes, batches, i, file_path, equilibrated, false);
  }


  update_cost();
  bool changed = false;
  it = 0;
  gap = 0.0;
  int inneriter = 1;
  double shift;
  int bush_index;

  // aon worker
  aonGraph dijfunc(this, dep, arr, dem, mode);

  // iterations
  if (verbose) Rcpp::Rcout << "Iterating..." <<  endl;

  while (it < max_iter){

    Rcpp::checkUserInterrupt();

    fill(equilibrated.begin(), equilibrated.end(), 0);

    for (int batch_index = 0; batch_index < n_batch; batch_index++){
      // load batch
      if (n_batch > 1) loadBushes(bushes, batches, batch_index, file_path, equilibrated, false);

      // optimize
      for (int i = 0; i <  batches[batch_index].size(); i++){


        Bush* b = &bushes[i];
        bush_index = batches[batch_index][i];

        b->root = od[bush_index].first;
        b->arr= od[bush_index].second.first;
        b->dem = od[bush_index].second.second;

        b->minmaxtree2(1);
        b->optimize2();
        if (b->changed) b->ordering();
        b->minmaxtree2(2);

        b->is_equilibrated= false;
        equilibrated[bush_index] = 0;

        double bushExcess = 0.0;

        for (int j = 0; j < nbedge; j++){

          if (bv.sdist[indG2[j]] < numeric_limits<double>::max() && bv.sdist[nodeG[j]] < numeric_limits<double>::max()) bushExcess += b->flow[j] * (wG[j] +bv.sdist[indG2[j]] - bv.sdist[nodeG[j]]);
        }

        double bushSPTT = 0.0;
        for (int j = 0; j < b->arr.size(); j++){
          if (bv.sdist[b->arr[j]] < numeric_limits<double>::max()) bushSPTT += bv.sdist[b->arr[j]] * b->dem[j];
        }


        //Rcpp::Rcout << bushExcess << " " <<  bushSPTT <<  " equilibrated!" << endl;
        if (bushExcess/bushSPTT < gap * 0.25) {
          b->is_equilibrated= true;
          equilibrated[bush_index] = 1;
          continue;
        }



        for (int j = nbnode - 1; j != -1; j--) {

          if (bv.slink[b->order[j]] == bv.llink[b->order[j]]) continue;
          if (bv.sparents[b->order[j]] == -1 ||  bv.lparents[b->order[j]] == -1) continue;
          if (abs(bv.sdist[b->order[j]] - bv.ldist[b->order[j]]) ==  0.0) continue;
          b->shift_flow(b->order[j]);

        }

      }
      // store batch
      if (n_batch > 1) storeBushes(bushes, batches, batch_index, file_path, equilibrated, false);
    }



      // equilibrate
      int count = 0;
      while (count < inneriter){

        for (int batch_index = 0; batch_index < n_batch; batch_index++){
          // load batch
          if (n_batch > 1) loadBushes(bushes, batches, batch_index, file_path, equilibrated, true);

        double count2 = 0;
        for (int i = 0; i < batches[batch_index].size(); i++){


          Bush* b = &bushes[i];
          bush_index = batches[batch_index][i];

          b->root = od[bush_index].first;
          b->arr= od[bush_index].second.first;
          b->dem = od[bush_index].second.second;

          //if (b->is_equilibrated) continue;
          if (equilibrated[bush_index] == 1) continue;

          b->minmaxtree2(2);

          double bushExcess = 0;
          for (int j = 0; j < nbedge; j++){

            if (bv.sdist[indG2[j]] < numeric_limits<double>::max() && bv.sdist[nodeG[j]] < numeric_limits<double>::max()) bushExcess += b->flow[j] * (wG[j] + bv.sdist[indG2[j]] - bv.sdist[nodeG[j]]);
          }

          double bushSPTT = 0.0;
          for (int j = 0; j < b->arr.size(); j++){
            if (bv.sdist[b->arr[j]] < numeric_limits<double>::max()) bushSPTT += bv.sdist[b->arr[j]] * b->dem[j];
          }


          if (bushExcess/bushSPTT < gap *0.25) count2 += 1;
          if (bushExcess/bushSPTT < gap * 0.25) {
            b->is_equilibrated = true;
            equilibrated[bush_index] = 1;
            continue;
          }

          for (int j = nbnode - 1; j != -1; j--) {

            if (bv.slink[b->order[j]] == bv.llink[b->order[j]]) continue;
            if (bv.sparents[b->order[j]] == -1 ||  bv.lparents[b->order[j]] == -1) continue;
            if (abs(bv.sdist[b->order[j]] - bv.ldist[b->order[j]]) == 0.0) continue;
            b->shift_flow(b->order[j]);

          }

        }



        // store batch
        if (n_batch > 1) storeBushes(bushes, batches, batch_index, file_path, equilibrated, true);
      }

        count += 1;
    }


    inneriter = inner_iter;

    it += 1;
    update_tstt();

    if (contract){
      aux = ch_aon(this, dep, arr, dem, phast, mode);
    }else{
      parallelReduce(0, dijfunc.n_iter, dijfunc);
      aux = dijfunc.m_result;
      fill(dijfunc.m_result.begin(), dijfunc.m_result.end(), 0.0);
    }


    update_sptt();
    gap = abs(tstt/sptt - 1);

    if (verbose) Rcpp::Rcout << "iteration " << it << " : " << gap <<  endl;
    if (gap < max_gap) break;


  }

}




