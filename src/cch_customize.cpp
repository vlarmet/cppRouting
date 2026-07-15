#include "cch_internal.h"

using namespace cch_internal;

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
                   IVec &backward_original){
  // Customization phase.
  //
  // The caller supplies one current cost per original edge. This function
  // updates the prepared shortcut graph so every CCH arc stores the best known
  // cost in each travel direction:
  // - forward: lower-ranked node -> higher-ranked node;
  // - backward: higher-ranked node -> lower-ranked node.
  //
  // For original arcs, *_original points to the input edge. For shortcuts,
  // *_first_arc and *_second_arc describe how to unpack the shortcut later.
  // CCH customization has two jobs:
  // 1. Load current edge weights onto the prepared CCH arcs.
  // 2. Relax all upward triangles so each shortcut stores the cheapest unpackable
  //    path for the current congestion state.
  forward.assign(tail.size(), INF);
  backward.assign(tail.size(), INF);
  forward_first_arc.assign(tail.size(), -1);
  forward_first_dir.assign(tail.size(), -1);
  forward_second_arc.assign(tail.size(), -1);
  forward_second_dir.assign(tail.size(), -1);
  forward_original.assign(tail.size(), -1);
  backward_first_arc.assign(tail.size(), -1);
  backward_first_dir.assign(tail.size(), -1);
  backward_second_arc.assign(tail.size(), -1);
  backward_second_dir.assign(tail.size(), -1);
  backward_original.assign(tail.size(), -1);

  for (int i = 0; i < gfrom.size(); i++){
    int arc = input_arc[i];
    if (arc < 0) continue;

    if (input_forward[i] == 1){
      if (gw[i] < forward[arc]){
        forward[arc] = gw[i];
        forward_original[arc] = i;
        forward_first_arc[arc] = -1;
        forward_second_arc[arc] = -1;
      }
    } else {
      if (gw[i] < backward[arc]){
        backward[arc] = gw[i];
        backward_original[arc] = i;
        backward_first_arc[arc] = -1;
        backward_second_arc[arc] = -1;
      }
    }
  }

  for (int lower_rank = 0; lower_rank < nb; lower_rank++){
    int begin = rank_first_out[lower_rank];
    int end = rank_first_out[lower_rank + 1];

    for (int i = begin; i < end; i++){
      int mid_rank = rank_adj_head[i];
      int left = rank_adj_arc[i];
      int mid_cursor = rank_first_out[mid_rank];

      for (int j = i + 1; j < end; j++){
        int upper_rank = rank_adj_head[j];
        int right = rank_adj_arc[j];
        int upper = find_rank_arc(mid_rank, upper_rank, rank_first_out, rank_adj_head, rank_adj_arc, mid_cursor);
        if (upper < 0) continue;

        // left and right form lower -> mid -> upper in rank space. Depending on
        // the original travel direction, this can improve either metric
        // direction stored on the upper shortcut.
        if (backward[left] != INF && forward[right] != INF){
          double candidate = backward[left] + forward[right];
          if (candidate < forward[upper]){
            forward[upper] = candidate;
            forward_original[upper] = -1;
            forward_first_arc[upper] = left;
            forward_first_dir[upper] = 0;
            forward_second_arc[upper] = right;
            forward_second_dir[upper] = 1;
          }
        }

        if (backward[right] != INF && forward[left] != INF){
          double candidate = backward[right] + forward[left];
          if (candidate < backward[upper]){
            backward[upper] = candidate;
            backward_original[upper] = -1;
            backward_first_arc[upper] = right;
            backward_first_dir[upper] = 0;
            backward_second_arc[upper] = left;
            backward_second_dir[upper] = 1;
          }
        }
      }
    }
  }
}
