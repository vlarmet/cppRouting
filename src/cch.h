#ifndef CCH_H
#define CCH_H

#include "graph.h"

struct CCHPrepared {
  int nbnode;

  // rank[node] is the contraction order position. Lower rank means the node is
  // contracted earlier. CCH arcs always point from lower rank to higher rank.
  IVec rank;

  // Topology of the prepared CCH, stored in original node ids. These arcs are
  // independent from edge weights and can be saved/reused while costs change.
  IVec tail;
  IVec head;
  IVec first_out;
  IVec adj_head;
  IVec adj_arc;

  // The same upward graph indexed by rank. Customization needs this form to
  // scan triangles in contraction order.
  IVec rank_first_out;
  IVec rank_adj_head;
  IVec rank_adj_arc;

  // Maps each input edge to the CCH arc that represents the same undirected
  // connection, plus the direction of the edge on that arc.
  IVec input_arc;
  IVec input_forward;

  // Parent pointers for the RoutingKit-style elimination-tree query. Walking
  // these parents visits the ancestors that can appear on an upward CCH search.
  IVec elimination_tree_parent;
};

CCHPrepared build_cch(IVec &gfrom, IVec &gto, int nb, IVec order);

void remap_cch_input_arcs(IVec &gfrom, IVec &gto, CCHPrepared &prepared);

void remap_cch_input_arcs(IVec &gfrom,
                          IVec &gto,
                          CCHPrepared &prepared,
                          IVec &input_arc,
                          IVec &input_forward);

void customize_cch(IVec &gfrom,
                   IVec &gto,
                   DVec &gw,
                   int nb,
                   IVec &rank,
                   IVec &tail,
                   IVec &head,
                   IVec &rank_first_out,
                   IVec &rank_adj_head,
                   IVec &rank_adj_arc,
                   IVec &input_arc,
                   IVec &input_forward,
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
                   IVec &backward_original);

DVec distance_pair_cch(int nb,
                       IVec &first_out,
                       IVec &adj_head,
                       IVec &adj_arc,
                       DVec &forward,
                       DVec &backward,
                       IVec dep,
                       IVec arr);

Rcpp::NumericMatrix distance_matrix_cch(int nb,
                                        IVec &first_out,
                                        IVec &adj_head,
                                        IVec &adj_arc,
                                        DVec &forward,
                                        DVec &backward,
                                        IVec dep,
                                        IVec arr);

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
                                         int n_values);

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
                                      DVec &demand);

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
                                       DVec demand);

#endif
