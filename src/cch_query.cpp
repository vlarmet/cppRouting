#include "cch_internal.h"

using namespace cch_internal;

struct CCHDistancePairWorker : public Worker {
  // One-to-one distance queries on a customized CCH metric. This is useful for
  // distance APIs and correctness checks, but it does not recover paths and is
  // not the traffic-assignment path.
  int nb;
  IVec &first_out;
  IVec &adj_head;
  IVec &adj_arc;
  DVec &forward;
  DVec &backward;
  IVec dep;
  IVec arr;
  DVec &result;

  CCHDistancePairWorker(int nb,
                        IVec &first_out,
                        IVec &adj_head,
                        IVec &adj_arc,
                        DVec &forward,
                        DVec &backward,
                        IVec dep,
                        IVec arr,
                        DVec &result)
    : nb(nb), first_out(first_out), adj_head(adj_head), adj_arc(adj_arc),
      forward(forward), backward(backward), dep(dep), arr(arr), result(result) {}

  void operator()(std::size_t begin, std::size_t end){
    DVec dist_forward(nb, INF);
    DVec dist_backward(nb, INF);
    IVec seen_forward(nb, 0);
    IVec seen_backward(nb, 0);
    IVec is_touched(nb, 0);
    IVec touched;
    touched.reserve(1024);

    for (std::size_t k = begin; k != end; k++){
      int source = dep[k];
      int target = arr[k];

      dist_forward[source] = 0.0;
      dist_backward[target] = 0.0;
      seen_forward[source] = 1;
      seen_backward[target] = 1;
      touch_node(source, is_touched, touched);
      touch_node(target, is_touched, touched);

      PQ forward_queue;
      PQ backward_queue;
      forward_queue.push(make_pair(source, 0.0));
      backward_queue.push(make_pair(target, 0.0));

      double best = source == target ? 0.0 : INF;

      while (!forward_queue.empty() || !backward_queue.empty()){
        double forward_min = forward_queue.empty() ? INF : forward_queue.top().second;
        double backward_min = backward_queue.empty() ? INF : backward_queue.top().second;
        if (forward_min >= best && backward_min >= best) break;

        if (!forward_queue.empty() && forward_min <= backward_min){
          int node = forward_queue.top().first;
          double dist = forward_queue.top().second;
          forward_queue.pop();

          if (dist > dist_forward[node]) continue;

          if (seen_backward[node] == 1){
            double candidate = dist_forward[node] + dist_backward[node];
            if (candidate < best) best = candidate;
          }

          for (int i = first_out[node]; i < first_out[node + 1]; i++){
            int arc = adj_arc[i];
            double edge_weight = forward[arc];
            if (edge_weight == INF) continue;

            int next = adj_head[i];
            double tentative = dist + edge_weight;
            if (tentative < dist_forward[next]){
              dist_forward[next] = tentative;
              seen_forward[next] = 1;
              touch_node(next, is_touched, touched);
              forward_queue.push(make_pair(next, tentative));

              if (seen_backward[next] == 1){
                double candidate = tentative + dist_backward[next];
                if (candidate < best) best = candidate;
              }
            }
          }
        } else if (!backward_queue.empty()){
          int node = backward_queue.top().first;
          double dist = backward_queue.top().second;
          backward_queue.pop();

          if (dist > dist_backward[node]) continue;

          if (seen_forward[node] == 1){
            double candidate = dist_forward[node] + dist_backward[node];
            if (candidate < best) best = candidate;
          }

          for (int i = first_out[node]; i < first_out[node + 1]; i++){
            int arc = adj_arc[i];
            double edge_weight = backward[arc];
            if (edge_weight == INF) continue;

            int next = adj_head[i];
            double tentative = dist + edge_weight;
            if (tentative < dist_backward[next]){
              dist_backward[next] = tentative;
              seen_backward[next] = 1;
              touch_node(next, is_touched, touched);
              backward_queue.push(make_pair(next, tentative));

              if (seen_forward[next] == 1){
                double candidate = dist_forward[next] + tentative;
                if (candidate < best) best = candidate;
              }
            }
          }
        }
      }

      result[k] = best;

      for (int i = 0; i < touched.size(); i++){
        int node = touched[i];
        dist_forward[node] = INF;
        dist_backward[node] = INF;
        seen_forward[node] = 0;
        seen_backward[node] = 0;
        is_touched[node] = 0;
      }
      touched.clear();
    }
  }
};

DVec distance_pair_cch(int nb,
                       IVec &first_out,
                       IVec &adj_head,
                       IVec &adj_arc,
                       DVec &forward,
                       DVec &backward,
                       IVec dep,
                       IVec arr){
  // Distance-only pair query. It computes shortest path costs but intentionally
  // does not recover original road edges, so it cannot be used for AON flow
  // loading.
  DVec result(dep.size(), INF);
  CCHDistancePairWorker worker(nb, first_out, adj_head, adj_arc, forward, backward, dep, arr, result);
  if (dep.size() <= 1000){
    worker(0, dep.size());
  } else {
    parallelFor(0, dep.size(), worker);
  }

  return result;
}

Rcpp::NumericMatrix distance_matrix_cch(int nb,
                                        IVec &first_out,
                                        IVec &adj_head,
                                        IVec &adj_arc,
                                        DVec &forward,
                                        DVec &backward,
                                        IVec dep,
                                        IVec arr){
  // Distance matrix helper. It precomputes one backward upward search per
  // destination, then combines it with each origin search.
  Rcpp::NumericMatrix result(dep.size(), arr.size());
  DVec dist_forward(nb, INF);
  DVec dist_backward(nb, INF);
  vector<DVec> backward_dist(arr.size(), DVec(nb, INF));

  for (int j = 0; j < arr.size(); j++){
    dijkstra_up(arr[j], nb, first_out, adj_head, adj_arc, backward, backward_dist[j]);
  }

  for (int i = 0; i < dep.size(); i++){
    dijkstra_up(dep[i], nb, first_out, adj_head, adj_arc, forward, dist_forward);
    for (int j = 0; j < arr.size(); j++){
      double best = INF;
      for (int node = 0; node < nb; node++){
        if (dist_forward[node] == INF || backward_dist[j][node] == INF) continue;
        double candidate = dist_forward[node] + backward_dist[j][node];
        if (candidate < best) best = candidate;
      }
      result(i, j) = best;
    }
    fill(dist_forward.begin(), dist_forward.end(), INF);
  }

  return result;
}

struct CCHPathValuesPairWorker : public Worker {
  // Pair query that routes on the customized CCH cost and accumulates one or
  // more original-edge value columns along the selected shortest path.
  IVec &gfrom;
  int nb;
  int n_values;
  IVec &rank;
  IVec &first_out;
  IVec &adj_head;
  IVec &adj_arc;
  IVec &elimination_tree_parent;
  DVec &forward;
  DVec &backward;
  IVec &forward_first_arc;
  IVec &forward_first_dir;
  IVec &forward_second_arc;
  IVec &forward_second_dir;
  IVec &forward_original;
  IVec &backward_first_arc;
  IVec &backward_first_dir;
  IVec &backward_second_arc;
  IVec &backward_second_dir;
  IVec &backward_original;
  IVec &dep;
  IVec &arr;
  DVec &values;
  vector<IVec> &groups;
  bool group_sources;
  DVec &result;

  CCHPathValuesPairWorker(IVec &gfrom,
                          int nb,
                          int n_values,
                          IVec &rank,
                          IVec &first_out,
                          IVec &adj_head,
                          IVec &adj_arc,
                          IVec &elimination_tree_parent,
                          DVec &forward,
                          DVec &backward,
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
                          IVec &dep,
                          IVec &arr,
                          DVec &values,
                          vector<IVec> &groups,
                          bool group_sources,
                          DVec &result)
    : gfrom(gfrom), nb(nb), n_values(n_values), rank(rank), first_out(first_out), adj_head(adj_head), adj_arc(adj_arc),
      elimination_tree_parent(elimination_tree_parent), forward(forward), backward(backward),
      forward_first_arc(forward_first_arc), forward_first_dir(forward_first_dir),
      forward_second_arc(forward_second_arc), forward_second_dir(forward_second_dir),
      forward_original(forward_original),
      backward_first_arc(backward_first_arc), backward_first_dir(backward_first_dir),
      backward_second_arc(backward_second_arc), backward_second_dir(backward_second_dir),
      backward_original(backward_original), dep(dep), arr(arr), values(values),
      groups(groups), group_sources(group_sources), result(result) {}

  CCHPathValuesPairWorker(CCHPathValuesPairWorker &other, Split)
    : gfrom(other.gfrom), nb(other.nb), n_values(other.n_values), rank(other.rank),
      first_out(other.first_out), adj_head(other.adj_head), adj_arc(other.adj_arc),
      elimination_tree_parent(other.elimination_tree_parent),
      forward(other.forward), backward(other.backward),
      forward_first_arc(other.forward_first_arc), forward_first_dir(other.forward_first_dir),
      forward_second_arc(other.forward_second_arc), forward_second_dir(other.forward_second_dir),
      forward_original(other.forward_original),
      backward_first_arc(other.backward_first_arc), backward_first_dir(other.backward_first_dir),
      backward_second_arc(other.backward_second_arc), backward_second_dir(other.backward_second_dir),
      backward_original(other.backward_original), dep(other.dep), arr(other.arr), values(other.values),
      groups(other.groups), group_sources(other.group_sources), result(other.result) {}

  void add_recovered_values(int pair_index,
                            int start,
                            int root,
                            bool reverse_path,
                            double cost,
                            IVec &pred_node,
                            IVec &pred_arc,
                            IVec &pred_dir,
                            IVec &path_arcs,
                            IVec &path_dirs,
                            IVec &path_edges){
    int n_pairs = dep.size();
    result[pair_index] = cost;
    for (int col = 0; col < n_values; col++) result[(col + 1) * n_pairs + pair_index] = 0.0;

    if (start == root) return;

    int node = start;
    while (node != root && pred_node[node] != -1){
      path_arcs.push_back(pred_arc[node]);
      path_dirs.push_back(pred_dir[node]);
      node = pred_node[node];
    }
    if (node != root){
      result[pair_index] = INF;
      for (int col = 0; col < n_values; col++) result[(col + 1) * n_pairs + pair_index] = INF;
      return;
    }
    if (reverse_path){
      reverse(path_arcs.begin(), path_arcs.end());
      reverse(path_dirs.begin(), path_dirs.end());
    }

    for (int i = 0; i < path_arcs.size(); i++){
      unpack_metric_arc(path_arcs[i], path_dirs[i],
                        forward_first_arc, forward_first_dir, forward_second_arc, forward_second_dir, forward_original,
                        backward_first_arc, backward_first_dir, backward_second_arc, backward_second_dir, backward_original,
                        path_edges);
    }
    for (int i = 0; i < path_edges.size(); i++){
      int edge = path_edges[i];
      if (edge < 0) continue;
      for (int col = 0; col < n_values; col++){
        result[(col + 1) * n_pairs + pair_index] += values[col * gfrom.size() + edge];
      }
    }
  }

  void operator()(std::size_t begin, std::size_t end){
    DVec dist(nb, INF);
    IVec pred_node(nb, -1);
    IVec pred_arc(nb, -1);
    IVec pred_dir(nb, -1);
    IVec is_touched(nb, 0);
    IVec touched;
    IVec is_pinned_ancestor(nb, 0);
    IVec pinned_stack;
    IVec pinned_touched;
    IVec path_arcs;
    IVec path_dirs;
    IVec path_edges;

    for (std::size_t group_id = begin; group_id != end; group_id++){
      if (groups[group_id].size() == 0) continue;

      int root = static_cast<int>(group_id);
      touch_distance_node(root, dist, pred_node, pred_arc, pred_dir, is_touched, touched);
      dist[root] = 0.0;

      if (group_sources){
        for_ancestors(elimination_tree_parent, root, [&](int node){
          relax_outgoing_elimination(node, first_out, adj_head, adj_arc, forward, 1,
                                     dist, pred_node, pred_arc, pred_dir, is_touched, touched);
          return true;
        });
      } else {
        for_ancestors(elimination_tree_parent, root, [&](int node){
          relax_outgoing_elimination(node, first_out, adj_head, adj_arc, backward, 0,
                                     dist, pred_node, pred_arc, pred_dir, is_touched, touched);
          return true;
        });
      }

      for (int pos = 0; pos < groups[group_id].size(); pos++){
        int k = groups[group_id][pos];
        int pinned = group_sources ? arr[k] : dep[k];
        for_ancestors(elimination_tree_parent, pinned, [&](int node){
          if (is_pinned_ancestor[node] == 0){
            is_pinned_ancestor[node] = 1;
            pinned_touched.push_back(node);
            pinned_stack.push_back(node);
          }
          return true;
        });
      }

      sort(pinned_stack.begin(), pinned_stack.end(), [&](int a, int b){
        return rank[a] > rank[b];
      });
      for (int stack_pos = 0; stack_pos < pinned_stack.size(); stack_pos++){
        int node = pinned_stack[stack_pos];
        if (group_sources){
          relax_incoming_elimination(node, first_out, adj_head, adj_arc, backward, 0,
                                     dist, pred_node, pred_arc, pred_dir, is_touched, touched);
        } else {
          relax_incoming_elimination(node, first_out, adj_head, adj_arc, forward, 1,
                                     dist, pred_node, pred_arc, pred_dir, is_touched, touched);
        }
      }

      for (int pos = 0; pos < groups[group_id].size(); pos++){
        int k = groups[group_id][pos];
        int start = group_sources ? arr[k] : dep[k];
        if (dist[start] == INF) continue;

        path_arcs.clear();
        path_dirs.clear();
        path_edges.clear();
        add_recovered_values(k, start, root, group_sources, dist[start],
                             pred_node, pred_arc, pred_dir, path_arcs, path_dirs, path_edges);
      }

      clear_elimination_state(dist, pred_node, pred_arc, pred_dir, is_touched, touched);
      for (int i = 0; i < pinned_touched.size(); i++) is_pinned_ancestor[pinned_touched[i]] = 0;
      pinned_touched.clear();
      pinned_stack.clear();
    }
  }
};

Rcpp::NumericMatrix path_values_pair_cch(IVec &gfrom,
                                         int nb,
                                         IVec &rank,
                                         IVec &first_out,
                                         IVec &adj_head,
                                         IVec &adj_arc,
                                         IVec &elimination_tree_parent,
                                         DVec &forward,
                                         DVec &backward,
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
                                         IVec dep,
                                         IVec arr,
                                         DVec &values,
                                         int n_values){
  IVec source_count(nb, 0);
  IVec target_count(nb, 0);
  int unique_sources = 0;
  int unique_targets = 0;
  for (int i = 0; i < dep.size(); i++){
    if (source_count[dep[i]] == 0) unique_sources++;
    if (target_count[arr[i]] == 0) unique_targets++;
    source_count[dep[i]]++;
    target_count[arr[i]]++;
  }

  bool group_sources = unique_sources <= unique_targets;
  vector<IVec> groups(nb);
  for (int i = 0; i < dep.size(); i++){
    int group_id = group_sources ? dep[i] : arr[i];
    groups[group_id].push_back(i);
  }

  DVec result((n_values + 1) * dep.size(), INF);
  CCHPathValuesPairWorker worker(gfrom, nb, n_values, rank, first_out, adj_head, adj_arc, elimination_tree_parent,
                                 forward, backward,
                                 forward_first_arc, forward_first_dir, forward_second_arc, forward_second_dir, forward_original,
                                 backward_first_arc, backward_first_dir, backward_second_arc, backward_second_dir, backward_original,
                                 dep, arr, values, groups, group_sources, result);
  if (dep.size() <= 1000 || nb <= 1000){
    worker(0, nb);
  } else {
    parallelFor(0, nb, worker);
  }

  Rcpp::NumericMatrix out(dep.size(), n_values + 1);
  for (int col = 0; col < n_values + 1; col++){
    for (int row = 0; row < dep.size(); row++){
      out(row, col) = result[col * dep.size() + row];
    }
  }
  return out;
}

struct CCHPathPairWorker : public Worker {
  // Pair query that recovers the original node path on a customized CCH metric.
  // It uses bidirectional Dijkstra on the CCH, then unpacks each CCH arc to
  // original edges and converts them to nodes via gfrom/gto.
  int nb;
  IVec &first_out;
  IVec &adj_head;
  IVec &adj_arc;
  DVec &forward;
  DVec &backward;
  IVec dep;
  IVec arr;
  IVec keep;
  vector<IVec> &result;
  // For unpacking CCH arcs to original edges and then to nodes
  IVec &gfrom;
  IVec &gto;
  IVec &forward_first_arc;
  IVec &forward_first_dir;
  IVec &forward_second_arc;
  IVec &forward_second_dir;
  IVec &forward_original;
  IVec &backward_first_arc;
  IVec &backward_first_dir;
  IVec &backward_second_arc;
  IVec &backward_second_dir;
  IVec &backward_original;

  CCHPathPairWorker(int nb,
                    IVec &first_out,
                    IVec &adj_head,
                    IVec &adj_arc,
                    DVec &forward,
                    DVec &backward,
                    IVec dep,
                    IVec arr,
                    IVec keep,
                    vector<IVec> &result,
                    IVec &gfrom,
                    IVec &gto,
                    IVec &forward_first_arc,
                    IVec &forward_first_dir,
                    IVec &forward_second_arc,
                    IVec &forward_second_dir,
                    IVec &forward_original,
                    IVec &backward_first_arc,
                    IVec &backward_first_dir,
                    IVec &backward_second_arc,
                    IVec &backward_second_dir,
                    IVec &backward_original)
    : nb(nb), first_out(first_out), adj_head(adj_head), adj_arc(adj_arc),
      forward(forward), backward(backward), dep(dep), arr(arr), keep(keep), result(result),
      gfrom(gfrom), gto(gto),
      forward_first_arc(forward_first_arc), forward_first_dir(forward_first_dir),
      forward_second_arc(forward_second_arc), forward_second_dir(forward_second_dir),
      forward_original(forward_original),
      backward_first_arc(backward_first_arc), backward_first_dir(backward_first_dir),
      backward_second_arc(backward_second_arc), backward_second_dir(backward_second_dir),
      backward_original(backward_original) {}

  CCHPathPairWorker(CCHPathPairWorker &other, Split)
    : nb(other.nb), first_out(other.first_out), adj_head(other.adj_head), adj_arc(other.adj_arc),
      forward(other.forward), backward(other.backward), dep(other.dep), arr(other.arr), keep(other.keep), result(other.result),
      gfrom(other.gfrom), gto(other.gto),
      forward_first_arc(other.forward_first_arc), forward_first_dir(other.forward_first_dir),
      forward_second_arc(other.forward_second_arc), forward_second_dir(other.forward_second_dir),
      forward_original(other.forward_original),
      backward_first_arc(other.backward_first_arc), backward_first_dir(other.backward_first_dir),
      backward_second_arc(other.backward_second_arc), backward_second_dir(other.backward_second_dir),
      backward_original(other.backward_original) {}

  void unpack_arcs_to_nodes(const IVec &path_arcs,
                            const IVec &path_dirs,
                            IVec &out_nodes){
    // Unpack CCH arcs to original edges, then convert to nodes.
    IVec path_edges;
    for (int i = 0; i < path_arcs.size(); i++){
      unpack_metric_arc(path_arcs[i], path_dirs[i],
                        forward_first_arc, forward_first_dir, forward_second_arc, forward_second_dir, forward_original,
                        backward_first_arc, backward_first_dir, backward_second_arc, backward_second_dir, backward_original,
                        path_edges);
    }

    // Convert original edges to nodes. For each edge (u,v), add v if it's not
    // the same as the last added node.
    for (int i = 0; i < path_edges.size(); i++){
      int edge = path_edges[i];
      if (edge < 0) continue;
      int u = gfrom[edge];
      int v = gto[edge];
      // We need to determine direction. Since unpack_metric_arc may have
      // reversed the edge, we check which endpoint matches the previous node.
      if (out_nodes.empty()){
        out_nodes.push_back(u);
        out_nodes.push_back(v);
      } else {
        int last = out_nodes.back();
        if (last == u){
          out_nodes.push_back(v);
        } else if (last == v){
          // Edge was traversed in reverse, add u
          out_nodes.push_back(u);
        }
      }
    }
  }

  void operator()(std::size_t begin, std::size_t end){
    DVec dist_forward(nb, INF);
    DVec dist_backward(nb, INF);
    IVec seen_forward(nb, 0);
    IVec seen_backward(nb, 0);
    IVec is_touched(nb, 0);
    IVec touched;
    IVec pred_forward(nb, -1);
    IVec pred_arc_forward(nb, -1);
    IVec pred_dir_forward(nb, -1);
    IVec pred_backward(nb, -1);
    IVec pred_arc_backward(nb, -1);
    IVec pred_dir_backward(nb, -1);

    for (std::size_t k = begin; k != end; k++){
      int start = dep[k];
      int end_node = arr[k];
      dist_forward[start] = 0.0;
      dist_backward[end_node] = 0.0;
      seen_forward[start] = 1;
      seen_backward[end_node] = 1;
      touch_node(start, is_touched, touched);
      touch_node(end_node, is_touched, touched);

      PQ forward_queue;
      PQ backward_queue;
      forward_queue.push(make_pair(start, 0.0));
      backward_queue.push(make_pair(end_node, 0.0));

      int mid = -1;
      double best = start == end_node ? 0.0 : INF;

      while (!forward_queue.empty() || !backward_queue.empty()){
        double forward_min = forward_queue.empty() ? INF : forward_queue.top().second;
        double backward_min = backward_queue.empty() ? INF : backward_queue.top().second;
        if (forward_min >= best && backward_min >= best) break;

        if (!forward_queue.empty() && forward_min <= backward_min){
          int node = forward_queue.top().first;
          double dist = forward_queue.top().second;
          forward_queue.pop();

          if (dist > dist_forward[node]) continue;

          if (seen_backward[node] == 1){
            double candidate = dist_forward[node] + dist_backward[node];
            if (candidate < best){
              best = candidate;
              mid = node;
            }
          }

          for (int i = first_out[node]; i < first_out[node + 1]; i++){
            int arc = adj_arc[i];
            double edge_weight = forward[arc];
            if (edge_weight == INF) continue;

            int next = adj_head[i];
            double tentative = dist + edge_weight;
            if (tentative < dist_forward[next]){
              dist_forward[next] = tentative;
              seen_forward[next] = 1;
              touch_node(next, is_touched, touched);
              forward_queue.push(make_pair(next, tentative));
              pred_forward[next] = node;
              pred_arc_forward[next] = arc;
              pred_dir_forward[next] = 1; // forward direction

              if (seen_backward[next] == 1){
                double candidate = tentative + dist_backward[next];
                if (candidate < best){
                  best = candidate;
                  mid = next;
                }
              }
            }
          }
        } else if (!backward_queue.empty()){
          int node = backward_queue.top().first;
          double dist = backward_queue.top().second;
          backward_queue.pop();

          if (dist > dist_backward[node]) continue;

          if (seen_forward[node] == 1){
            double candidate = dist_forward[node] + dist_backward[node];
            if (candidate < best){
              best = candidate;
              mid = node;
            }
          }

          for (int i = first_out[node]; i < first_out[node + 1]; i++){
            int arc = adj_arc[i];
            double edge_weight = backward[arc];
            if (edge_weight == INF) continue;

            int next = adj_head[i];
            double tentative = dist + edge_weight;
            if (tentative < dist_backward[next]){
              dist_backward[next] = tentative;
              seen_backward[next] = 1;
              touch_node(next, is_touched, touched);
              backward_queue.push(make_pair(next, tentative));
              pred_backward[next] = node;
              pred_arc_backward[next] = arc;
              pred_dir_backward[next] = 0; // backward direction

              if (seen_forward[next] == 1){
                double candidate = dist_forward[next] + tentative;
                if (candidate < best){
                  best = candidate;
                  mid = next;
                }
              }
            }
          }
        }
      }

      IVec path_result;
      if (mid != -1){
        // Reconstruct path as sequence of CCH arcs with directions using pred_arc/pred_dir
        IVec path_arcs;
        IVec path_dirs;

        // Forward part: from start to mid (follow pred_forward backwards)
        {
          int node = mid;
          while (node != start && pred_forward[node] != -1){
            path_arcs.push_back(pred_arc_forward[node]);
            path_dirs.push_back(pred_dir_forward[node]);
            node = pred_forward[node];
          }
          // Reverse to get forward order
          reverse(path_arcs.begin(), path_arcs.end());
          reverse(path_dirs.begin(), path_dirs.end());
        }

        // Backward part: from mid to end_node (follow pred_backward forwards)
        {
          int node = mid;
          while (node != end_node && pred_backward[node] != -1){
            path_arcs.push_back(pred_arc_backward[node]);
            path_dirs.push_back(pred_dir_backward[node]);
            node = pred_backward[node];
          }
        }

        // Unpack CCH arcs to original edges, then to nodes
        unpack_arcs_to_nodes(path_arcs, path_dirs, path_result);
      }

      result[k] = path_result;
    }
  }
};

vector<IVec> path_pair_cch(int nb,
                           IVec &first_out,
                           IVec &adj_head,
                           IVec &adj_arc,
                           DVec &forward,
                           DVec &backward,
                           IVec &gfrom,
                           IVec &gto,
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
                           IVec dep,
                           IVec arr,
                           IVec keep){
  vector<IVec> result(dep.size());
  CCHPathPairWorker worker(nb, first_out, adj_head, adj_arc, forward, backward, dep, arr, keep, result,
                           gfrom, gto, forward_first_arc, forward_first_dir, forward_second_arc, forward_second_dir, forward_original,
                           backward_first_arc, backward_first_dir, backward_second_arc, backward_second_dir, backward_original);
  if (dep.size() <= 1000){
    worker(0, dep.size());
  } else {
    parallelFor(0, dep.size(), worker);
  }
  return result;
}
