#include "cch_internal.h"

using namespace cch_internal;

CCHPrepared build_cch(IVec &gfrom, IVec &gto, int nb, IVec order){
  // Preparation phase.
  //
  // Input:
  // - original directed road edges, using integer node ids;
  // - an optional node order.
  //
  // Output:
  // - a directed upward shortcut graph;
  // - the elimination tree used by fast grouped queries;
  // - a map from each original road edge to the CCH arc where its cost enters.
  //
  // No travel-time/congestion weight is stored here. That is intentional: the
  // result can be saved and reused as long as the road topology is unchanged.
  if (order.size() == 0) order = min_degree_order(gfrom, gto, nb);

  CCHPrepared result;
  result.nbnode = nb;
  result.rank = rank_from_order(order, nb);

  // Build the upward graph induced by the order. Each input edge is first
  // treated as undirected; direction is restored during customization.
  vector<IVec> upward(nb);
  for (int i = 0; i < gfrom.size(); i++){
    int from_rank = result.rank[gfrom[i]];
    int to_rank = result.rank[gto[i]];
    if (from_rank == to_rank) continue;
    if (from_rank < to_rank){
      upward[from_rank].push_back(to_rank);
    } else {
      upward[to_rank].push_back(from_rank);
    }
  }

  for (int node_rank = 0; node_rank < nb; node_rank++) sort_unique(upward[node_rank]);

  vector<pair<int, int> > rank_arcs;
  rank_arcs.reserve(gfrom.size());

  for (int node_rank = 0; node_rank < nb; node_rank++){
    IVec &neighbors = upward[node_rank];
    if (neighbors.size() == 0) continue;

    for (int i = 0; i < neighbors.size(); i++) rank_arcs.push_back(make_pair(node_rank, neighbors[i]));

    // Eliminate this rank by making all higher-ranked neighbors adjacent. The
    // fill edges are what let later customization evaluate shortcut triangles.
    int lowest_upward_neighbor = neighbors[0];
    if (neighbors.size() > 1){
      IVec fill_edges(neighbors.begin() + 1, neighbors.end());
      upward[lowest_upward_neighbor] = merge_unique_tail(upward[lowest_upward_neighbor], fill_edges);
    }
  }

  sort(rank_arcs.begin(), rank_arcs.end());
  rank_arcs.erase(unique(rank_arcs.begin(), rank_arcs.end()), rank_arcs.end());

  result.tail.reserve(rank_arcs.size());
  result.head.reserve(rank_arcs.size());

  vector<IVec> cch_upward(nb);
  for (int arc = 0; arc < rank_arcs.size(); arc++){
    int lower_rank = rank_arcs[arc].first;
    int higher_rank = rank_arcs[arc].second;
    int lower_node = order[lower_rank];
    int higher_node = order[higher_rank];

    result.tail.push_back(lower_node);
    result.head.push_back(higher_node);
    cch_upward[lower_rank].push_back(higher_rank);
  }

  result.rank_first_out.assign(nb + 1, 0);
  for (int rank_id = 0; rank_id < nb; rank_id++) result.rank_first_out[rank_id + 1] = cch_upward[rank_id].size();
  for (int rank_id = 1; rank_id <= nb; rank_id++) result.rank_first_out[rank_id] += result.rank_first_out[rank_id - 1];

  result.rank_adj_head.resize(result.tail.size());
  result.rank_adj_arc.resize(result.tail.size());
  IVec rank_cursor = result.rank_first_out;
  for (int arc = 0; arc < rank_arcs.size(); arc++){
    int lower_rank = rank_arcs[arc].first;
    int higher_rank = rank_arcs[arc].second;
    int pos = rank_cursor[lower_rank]++;
    result.rank_adj_head[pos] = higher_rank;
    result.rank_adj_arc[pos] = arc;
  }

  result.first_out.assign(nb + 1, 0);
  for (int i = 0; i < result.tail.size(); i++) result.first_out[result.tail[i] + 1]++;
  for (int i = 1; i <= nb; i++) result.first_out[i] += result.first_out[i - 1];

  result.adj_head.resize(result.tail.size());
  result.adj_arc.resize(result.tail.size());
  IVec cursor = result.first_out;
  for (int arc = 0; arc < result.tail.size(); arc++){
    int pos = cursor[result.tail[arc]]++;
    result.adj_head[pos] = result.head[arc];
    result.adj_arc[pos] = arc;
  }
  result.elimination_tree_parent = elimination_tree_from_up_graph(nb, result.first_out, result.adj_head);

  // Remember where each original edge enters the CCH. Congestion assignment
  // changes edge costs many times, so customization must not search for this
  // mapping at every iteration.
  unordered_map<long long, int> edge_to_arc;
  edge_to_arc.reserve(result.tail.size() * 2);
  for (int arc = 0; arc < result.tail.size(); arc++){
    edge_to_arc[edge_key(result.tail[arc], result.head[arc])] = arc;
  }

  result.input_arc.assign(gfrom.size(), -1);
  result.input_forward.assign(gfrom.size(), 0);
  for (int i = 0; i < gfrom.size(); i++){
    int from = gfrom[i];
    int to = gto[i];
    if (from == to) continue;

    int lower = result.rank[from] < result.rank[to] ? from : to;
    int higher = result.rank[from] < result.rank[to] ? to : from;
    unordered_map<long long, int>::iterator found = edge_to_arc.find(edge_key(lower, higher));
    if (found == edge_to_arc.end()) continue;

    result.input_arc[i] = found->second;
    result.input_forward[i] = (from == lower && to == higher) ? 1 : 0;
  }

  return result;
}

void remap_cch_input_arcs(IVec &gfrom, IVec &gto, CCHPrepared &prepared){
  remap_cch_input_arcs(gfrom, gto, prepared, prepared.input_arc, prepared.input_forward);
}

void remap_cch_input_arcs(IVec &gfrom,
                          IVec &gto,
                          CCHPrepared &prepared,
                          IVec &input_arc,
                          IVec &input_forward){
  // Prepared CCH objects can be reused by assign_traffic(). The assignment
  // graph has the same topology but its current adjacency edge order may differ,
  // so only this cheap input-edge map is rebuilt.
  unordered_map<long long, int> edge_to_arc;
  edge_to_arc.reserve(prepared.tail.size() * 2);
  for (int arc = 0; arc < prepared.tail.size(); arc++){
    edge_to_arc[edge_key(prepared.tail[arc], prepared.head[arc])] = arc;
  }

  input_arc.assign(gfrom.size(), -1);
  input_forward.assign(gfrom.size(), 0);
  for (int i = 0; i < gfrom.size(); i++){
    int from = gfrom[i];
    int to = gto[i];
    if (from == to) continue;

    int lower = prepared.rank[from] < prepared.rank[to] ? from : to;
    int higher = prepared.rank[from] < prepared.rank[to] ? to : from;
    unordered_map<long long, int>::iterator found = edge_to_arc.find(edge_key(lower, higher));
    if (found == edge_to_arc.end()) continue;

    input_arc[i] = found->second;
    input_forward[i] = (from == lower && to == higher) ? 1 : 0;
  }
}
