#ifndef CCH_INTERNAL_H
#define CCH_INTERNAL_H

#include "cch.h"
#include <algorithm>
#include <limits>
#include <queue>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <progress.hpp>
#include <progress_bar.hpp>

using namespace std;
using namespace RcppParallel;

namespace cch_internal {

const double INF = numeric_limits<double>::max();

inline long long edge_key(int a, int b){
  if (a > b) std::swap(a, b);
  return (static_cast<long long>(a) << 32) | static_cast<unsigned int>(b);
}

inline IVec identity_order(int nb){
  IVec order(nb);
  for (int i = 0; i < nb; i++) order[i] = i;
  return order;
}

// [[Rcpp::depends(RcppProgress)]]
inline IVec min_degree_order(IVec &gfrom, IVec &gto, int nb, bool display_progress = false){
  vector<set<int> > graph(nb);
  for (int i = 0; i < gfrom.size(); i++){
    if (gfrom[i] == gto[i]) continue;
    graph[gfrom[i]].insert(gto[i]);
    graph[gto[i]].insert(gfrom[i]);
  }

  IVec order;
  order.reserve(nb);
  IVec removed(nb, 0);

  // Bucket-by-degree queue: buckets[d] holds, in ascending node-id order, the
  // unremoved nodes whose current degree is d. This finds the minimum-degree
  // node (ties broken by smallest id, same as a linear scan) in amortized
  // near-O(1) instead of rescanning every node at each of the nb elimination
  // steps, which made preparation O(nb^2) and dominated cpp_cch_prepare's
  // runtime on large road networks.
  IVec degree(nb);
  int max_degree = 0;
  for (int node = 0; node < nb; node++){
    degree[node] = (int)graph[node].size();
    max_degree = max(max_degree, degree[node]);
  }

  vector<set<int> > buckets(max_degree + 1);
  for (int node = 0; node < nb; node++) buckets[degree[node]].insert(node);

  int min_ptr = 0;

  Progress p(nb, display_progress);
  for (int step = 0; step < nb; step++){
    Rcpp::checkUserInterrupt();
    p.increment();

    while (min_ptr < (int)buckets.size() && buckets[min_ptr].empty()) min_ptr++;
    if (min_ptr >= (int)buckets.size()) break;

    int best = *buckets[min_ptr].begin();
    buckets[min_ptr].erase(buckets[min_ptr].begin());

    IVec neighbors;
    for (set<int>::iterator it = graph[best].begin(); it != graph[best].end(); ++it){
      if (removed[*it] == 0) neighbors.push_back(*it);
    }

    for (int i = 0; i < neighbors.size(); i++){
      for (int j = i + 1; j < neighbors.size(); j++){
        graph[neighbors[i]].insert(neighbors[j]);
        graph[neighbors[j]].insert(neighbors[i]);
      }
    }

    for (int i = 0; i < neighbors.size(); i++) graph[neighbors[i]].erase(best);
    graph[best].clear();
    removed[best] = 1;
    order.push_back(best);

    // Fill-in only changes the degree of best's former neighbors: they lose
    // the arc to best and gain any new clique edges among themselves.
    for (int i = 0; i < neighbors.size(); i++){
      int node = neighbors[i];
      int new_degree = (int)graph[node].size();
      if (new_degree == degree[node]) continue;

      buckets[degree[node]].erase(node);
      if (new_degree >= (int)buckets.size()) buckets.resize(new_degree + 1);
      buckets[new_degree].insert(node);
      degree[node] = new_degree;
      if (new_degree < min_ptr) min_ptr = new_degree;
    }
  }

  return order;
}

inline IVec rank_from_order(IVec order, int nb){
  if (order.size() == 0) order = identity_order(nb);
  if (order.size() != nb) throw invalid_argument("order length must be equal to the number of nodes");

  IVec rank(nb, -1);
  for (int i = 0; i < nb; i++){
    int node = order[i];
    if (node < 0 || node >= nb) throw invalid_argument("order contains an invalid node id");
    if (rank[node] != -1) throw invalid_argument("order contains duplicated node ids");
    rank[node] = i;
  }
  return rank;
}

inline void dijkstra_up(int start,
                 int nb,
                 IVec &first_out,
                 IVec &adj_head,
                 IVec &adj_arc,
                 DVec &weight,
                 DVec &distances,
                 IVec &touched){
  // Unidirectional Dijkstra restricted to the upward CCH graph. Also records
  // which nodes were actually reached, so callers can combine many such
  // searches (e.g. a many-to-many distance matrix) without rescanning the
  // full node set for every pair.
  distances[start] = 0.0;
  touched.push_back(start);
  PQ queue;
  queue.push(make_pair(start, 0.0));

  while (!queue.empty()){
    int node = queue.top().first;
    double dist = queue.top().second;
    queue.pop();

    if (dist > distances[node]) continue;

    for (int i = first_out[node]; i < first_out[node + 1]; i++){
      int arc = adj_arc[i];
      int next = adj_head[i];
      double edge_weight = weight[arc];
      if (edge_weight == INF) continue;

      double tentative = dist + edge_weight;
      if (tentative < distances[next]){
        if (distances[next] == INF) touched.push_back(next);
        distances[next] = tentative;
        queue.push(make_pair(next, tentative));
      }
    }
  }
}

inline void touch_node(int node, IVec &is_touched, IVec &touched){
  if (is_touched[node] == 0){
    is_touched[node] = 1;
    touched.push_back(node);
  }
}

inline void sort_unique(IVec &values){
  sort(values.begin(), values.end());
  values.erase(unique(values.begin(), values.end()), values.end());
}

inline IVec merge_unique_tail(IVec &base, IVec &added){
  IVec result;
  result.reserve(base.size() + added.size());
  merge(base.begin(), base.end(), added.begin(), added.end(), back_inserter(result));
  result.erase(unique(result.begin(), result.end()), result.end());
  return result;
}

inline int find_rank_arc(int tail_rank,
                  int head_rank,
                  IVec &rank_first_out,
                  IVec &rank_adj_head,
                  IVec &rank_adj_arc,
                  int &cursor){
  int end = rank_first_out[tail_rank + 1];
  while (cursor < end && rank_adj_head[cursor] < head_rank) cursor++;
  if (cursor < end && rank_adj_head[cursor] == head_rank) return rank_adj_arc[cursor];
  return -1;
}

inline IVec elimination_tree_from_up_graph(int nb, IVec &first_out, IVec &adj_head){
  IVec parent(nb, -1);
  for (int node = 0; node < nb; node++){
    if (first_out[node] != first_out[node + 1]) parent[node] = adj_head[first_out[node]];
  }
  return parent;
}

template<class F>
inline void for_ancestors(IVec &parent, int node, const F &f){
  while (node != -1){
    if (!f(node)) return;
    node = parent[node];
  }
}

inline void touch_distance_node(int node, DVec &dist, IVec &pred_node, IVec &pred_arc, IVec &pred_dir,
                         IVec &is_touched, IVec &touched){
  if (is_touched[node] == 0){
    is_touched[node] = 1;
    touched.push_back(node);
  }
}

inline void relax_outgoing_elimination(int node,
                                IVec &first_out,
                                IVec &adj_head,
                                IVec &adj_arc,
                                DVec &weight,
                                int arc_dir,
                                DVec &dist,
                                IVec &pred_node,
                                IVec &pred_arc,
                                IVec &pred_dir,
                                IVec &is_touched,
                                IVec &touched){
  if (dist[node] == INF) return;
  for (int i = first_out[node]; i < first_out[node + 1]; i++){
    int arc = adj_arc[i];
    double edge_weight = weight[arc];
    if (edge_weight == INF) continue;

    int next = adj_head[i];
    double tentative = dist[node] + edge_weight;
    if (tentative < dist[next]){
      touch_distance_node(next, dist, pred_node, pred_arc, pred_dir, is_touched, touched);
      dist[next] = tentative;
      pred_node[next] = node;
      pred_arc[next] = arc;
      pred_dir[next] = arc_dir;
    }
  }
}

inline void relax_incoming_elimination(int node,
                                IVec &first_out,
                                IVec &adj_head,
                                IVec &adj_arc,
                                DVec &weight,
                                int arc_dir,
                                DVec &dist,
                                IVec &pred_node,
                                IVec &pred_arc,
                                IVec &pred_dir,
                                IVec &is_touched,
                                IVec &touched){
  for (int i = first_out[node]; i < first_out[node + 1]; i++){
    int arc = adj_arc[i];
    double edge_weight = weight[arc];
    if (edge_weight == INF) continue;

    int next = adj_head[i];
    if (dist[next] == INF) continue;
    double tentative = dist[next] + edge_weight;
    if (tentative < dist[node]){
      touch_distance_node(node, dist, pred_node, pred_arc, pred_dir, is_touched, touched);
      dist[node] = tentative;
      pred_node[node] = next;
      pred_arc[node] = arc;
      pred_dir[node] = arc_dir;
    }
  }
}

inline void clear_elimination_state(DVec &dist,
                             IVec &pred_node,
                             IVec &pred_arc,
                             IVec &pred_dir,
                             IVec &is_touched,
                             IVec &touched){
  for (int i = 0; i < touched.size(); i++){
    int node = touched[i];
    dist[node] = INF;
    pred_node[node] = -1;
    pred_arc[node] = -1;
    pred_dir[node] = -1;
    is_touched[node] = 0;
  }
  touched.clear();
}

inline void unpack_metric_arc(int arc,
                       int direction,
                       IVec &forward_first_arc,
                       IVec &forward_first_dir,
                       IVec &forward_second_arc,
                       IVec &forward_second_dir,
                       IVec &forward_original,
                       IVec &backward_first_arc,
                       IVec &backward_first_dir,
                       IVec &backward_second_arc,
                       IVec &backward_second_dir,
                       IVec &backward_original,
                       IVec &path_edges){
  if (arc < 0) return;

  if (direction == 1){
    if (forward_original[arc] >= 0){
      path_edges.push_back(forward_original[arc]);
      return;
    }
    unpack_metric_arc(forward_first_arc[arc], forward_first_dir[arc],
                      forward_first_arc, forward_first_dir, forward_second_arc, forward_second_dir, forward_original,
                      backward_first_arc, backward_first_dir, backward_second_arc, backward_second_dir, backward_original,
                      path_edges);
    unpack_metric_arc(forward_second_arc[arc], forward_second_dir[arc],
                      forward_first_arc, forward_first_dir, forward_second_arc, forward_second_dir, forward_original,
                      backward_first_arc, backward_first_dir, backward_second_arc, backward_second_dir, backward_original,
                      path_edges);
  } else {
    if (backward_original[arc] >= 0){
      path_edges.push_back(backward_original[arc]);
      return;
    }
    unpack_metric_arc(backward_first_arc[arc], backward_first_dir[arc],
                      forward_first_arc, forward_first_dir, forward_second_arc, forward_second_dir, forward_original,
                      backward_first_arc, backward_first_dir, backward_second_arc, backward_second_dir, backward_original,
                      path_edges);
    unpack_metric_arc(backward_second_arc[arc], backward_second_dir[arc],
                      forward_first_arc, forward_first_dir, forward_second_arc, forward_second_dir, forward_original,
                      backward_first_arc, backward_first_dir, backward_second_arc, backward_second_dir, backward_original,
                      path_edges);
  }
}

} // namespace cch_internal

#endif
