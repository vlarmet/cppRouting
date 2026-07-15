#include "cch_internal.h"

using namespace cch_internal;

struct CCHEliminationGroupedAonWorker : public Worker {
  // RoutingKit-style grouped elimination-tree AON query.
  //
  // Instead of running one independent shortest-path search per OD pair, this
  // worker groups OD pairs by the side with fewer unique endpoints. A group then
  // shares one search from that repeated endpoint and only performs the smaller
  // pinned-endpoint work needed for the individual OD pairs.
  //
  // The recovered paths are unpacked to original graph edges immediately and
  // accumulated in flow. That is why this worker is the traffic-assignment CCH
  // path, while the distance query above is not.
  IVec &gfrom;
  int nb;
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
  DVec &demand;
  vector<IVec> &groups;
  bool group_sources;
  DVec flow;

  CCHEliminationGroupedAonWorker(IVec &gfrom,
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
                                 IVec &dep,
                                 IVec &arr,
                                 DVec &demand,
                                 vector<IVec> &groups,
                                 bool group_sources)
    : gfrom(gfrom), nb(nb), rank(rank), first_out(first_out), adj_head(adj_head), adj_arc(adj_arc),
      elimination_tree_parent(elimination_tree_parent), forward(forward), backward(backward),
      forward_first_arc(forward_first_arc), forward_first_dir(forward_first_dir),
      forward_second_arc(forward_second_arc), forward_second_dir(forward_second_dir),
      forward_original(forward_original),
      backward_first_arc(backward_first_arc), backward_first_dir(backward_first_dir),
      backward_second_arc(backward_second_arc), backward_second_dir(backward_second_dir),
      backward_original(backward_original), dep(dep), arr(arr), demand(demand),
      groups(groups), group_sources(group_sources), flow(gfrom.size(), 0.0) {}

  CCHEliminationGroupedAonWorker(CCHEliminationGroupedAonWorker &other, Split)
    : gfrom(other.gfrom), nb(other.nb), rank(other.rank),
      first_out(other.first_out), adj_head(other.adj_head), adj_arc(other.adj_arc),
      elimination_tree_parent(other.elimination_tree_parent),
      forward(other.forward), backward(other.backward),
      forward_first_arc(other.forward_first_arc), forward_first_dir(other.forward_first_dir),
      forward_second_arc(other.forward_second_arc), forward_second_dir(other.forward_second_dir),
      forward_original(other.forward_original),
      backward_first_arc(other.backward_first_arc), backward_first_dir(other.backward_first_dir),
      backward_second_arc(other.backward_second_arc), backward_second_dir(other.backward_second_dir),
      backward_original(other.backward_original), dep(other.dep), arr(other.arr), demand(other.demand),
      groups(other.groups), group_sources(other.group_sources), flow(other.gfrom.size(), 0.0) {}

  void join(CCHEliminationGroupedAonWorker &other){
    for (int i = 0; i < flow.size(); i++) flow[i] += other.flow[i];
  }

  void add_recovered_path(int start,
                          int root,
                          bool reverse_path,
                          double volume,
                          IVec &pred_node,
                          IVec &pred_arc,
                          IVec &pred_dir,
                          IVec &path_arcs,
                          IVec &path_dirs,
                          IVec &path_edges){
    int node = start;
    while (node != root && pred_node[node] != -1){
      path_arcs.push_back(pred_arc[node]);
      path_dirs.push_back(pred_dir[node]);
      node = pred_node[node];
    }
    if (node != root) {
      path_arcs.clear();
      path_dirs.clear();
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
      if (path_edges[i] >= 0) flow[path_edges[i]] += volume;
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

      // Shared phase: one upward search from the repeated endpoint. For source
      // groups this is a forward search from the origin; for destination groups
      // it is the symmetric backward search from the destination.
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

      // Pinned phase: only visit ancestors of the OD endpoints in this group.
      // Processing high ranks first is essential. It makes the dynamic program
      // propagate distances from common ancestors back down toward each pinned
      // endpoint before we recover paths.
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

      // Recover and unpack one shortest path for every OD pair in the group.
      // Shortcut unpacking writes flow back to original input edges, which is
      // what traffic assignment needs for the next cost update.
      for (int pos = 0; pos < groups[group_id].size(); pos++){
        int k = groups[group_id][pos];
        int start = group_sources ? arr[k] : dep[k];
        if (start == root || dist[start] == INF) continue;

        path_arcs.clear();
        path_dirs.clear();
        path_edges.clear();
        add_recovered_path(start, root, group_sources, demand[k],
                           pred_node, pred_arc, pred_dir, path_arcs, path_dirs, path_edges);
      }

      clear_elimination_state(dist, pred_node, pred_arc, pred_dir, is_touched, touched);
      for (int i = 0; i < pinned_touched.size(); i++) is_pinned_ancestor[pinned_touched[i]] = 0;
      pinned_touched.clear();
      pinned_stack.clear();
    }
  }
};

DVec aon_flow_cch_elimination_grouped(IVec &gfrom,
                                      IVec &gto,
                                      DVec &gw,
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
                                      IVec &dep,
                                      IVec &arr,
                                      DVec &demand){
  // Public native AON entry point for the current CCH implementation. The R API
  // exposes this as aon_method = "cch". Older pairwise/grouped CCH query modes
  // were removed so assignment always uses this elimination-tree path.
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

  CCHEliminationGroupedAonWorker worker(gfrom, nb, rank, first_out, adj_head, adj_arc, elimination_tree_parent,
                                        forward, backward,
                                        forward_first_arc, forward_first_dir, forward_second_arc, forward_second_dir, forward_original,
                                        backward_first_arc, backward_first_dir, backward_second_arc, backward_second_dir, backward_original,
                                        dep, arr, demand, groups, group_sources);
  if (dep.size() <= 1000 || nb <= 1000){
    worker(0, nb);
  } else {
    parallelReduce(0, nb, worker);
  }

  return worker.flow;
}

Rcpp::List aon_cch_elimination_grouped(IVec &gfrom,
                                       IVec &gto,
                                       DVec &gw,
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
                                       DVec demand){
  // Thin Rcpp-friendly wrapper around the native flow vector. It returns the
  // same four-column shape as the other cppRouting AON backends.
  DVec flow = aon_flow_cch_elimination_grouped(gfrom, gto, gw, nb, rank, first_out, adj_head, adj_arc,
                                               elimination_tree_parent, forward, backward,
                                               forward_first_arc, forward_first_dir, forward_second_arc, forward_second_dir, forward_original,
                                               backward_first_arc, backward_first_dir, backward_second_arc, backward_second_dir, backward_original,
                                               dep, arr, demand);

  Rcpp::List result(4);
  result[0] = gfrom;
  result[1] = gto;
  result[2] = gw;
  result[3] = flow;
  return result;
}
