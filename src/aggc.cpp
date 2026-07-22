#include <vector>
#include <RcppParallel.h>
#include <Rcpp.h>
#include "graph.h"
#include "cgraph.h"
#include "aggc.h"
using namespace RcppParallel;

aggC::aggC(CGraph* gr, Graph* original):
  m_gr(gr), m_or(original){

  m_result.resize(m_gr->nodeG.size(), 0.0);
  m_result2.resize(m_gr->nodeGr.size(), 0.0);
}

void aggC::operator()(std::size_t begin, std::size_t end){

  for (size_t k = begin; k != end; k++){
    // upgraph
    for (int i = m_gr->indG[k]; i < m_gr->indG[k+1]; i++){

      double weight = 0.0;
      IVec result2(2);
      result2[0] = k;
      result2[1] = m_gr->nodeG[i];
      m_gr->unpack(result2);

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

        if (edge_index == -1) Rcpp::stop("aggC: no matching original edge found for %d->%d", dep, second);

        weight += m_or->add[edge_index];
        dep = second;
        pos += 1;
      }

        m_result[i] = weight;
    }

    // downgraph
    for (int i = m_gr->indGr[k]; i < m_gr->indGr[k+1]; i++){

      double weight = 0.0;
      IVec result2(2);
      result2[0] = m_gr->nodeGr[i];;
      result2[1] = k;
      m_gr->unpack(result2);

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

        if (edge_index == -1) Rcpp::stop("aggC: no matching original edge found for %d->%d", dep, second);

        weight += m_or->add[edge_index];
        dep = second;
        pos += 1;
      }

      m_result2[i] = weight;
    }

  }
}

aggCAdd::aggCAdd(CGraph* gr, Graph* original):
  m_gr(gr), m_or(original){

  n_add = m_or->n_add;
  m_result.assign((size_t)m_gr->nodeG.size() * n_add, 0.0);
  m_result2.assign((size_t)m_gr->nodeGr.size() * n_add, 0.0);
}

void aggCAdd::operator()(std::size_t begin, std::size_t end){

  int total = (int)m_or->nodeG.size(); // stride of m_or->add
  int strideup = (int)m_gr->nodeG.size();
  int stridedown = (int)m_gr->nodeGr.size();

  for (size_t k = begin; k != end; k++){
    // upgraph
    for (int i = m_gr->indG[k]; i < m_gr->indG[k+1]; i++){

      DVec weight(n_add, 0.0);
      IVec result2(2);
      result2[0] = k;
      result2[1] = m_gr->nodeG[i];
      m_gr->unpack(result2);

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

        if (edge_index == -1) Rcpp::stop("aggCAdd: no matching original edge found for %d->%d", dep, second);

        for (int col = 0; col < n_add; col++){
          weight[col] += m_or->add[(size_t)col * total + edge_index];
        }
        dep = second;
        pos += 1;
      }

      for (int col = 0; col < n_add; col++){
        m_result[(size_t)col * strideup + i] = weight[col];
      }
    }

    // downgraph
    for (int i = m_gr->indGr[k]; i < m_gr->indGr[k+1]; i++){

      DVec weight(n_add, 0.0);
      IVec result2(2);
      result2[0] = m_gr->nodeGr[i];;
      result2[1] = k;
      m_gr->unpack(result2);

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

        if (edge_index == -1) Rcpp::stop("aggCAdd: no matching original edge found for %d->%d", dep, second);

        for (int col = 0; col < n_add; col++){
          weight[col] += m_or->add[(size_t)col * total + edge_index];
        }
        dep = second;
        pos += 1;
      }

      for (int col = 0; col < n_add; col++){
        m_result2[(size_t)col * stridedown + i] = weight[col];
      }
    }

  }
}
