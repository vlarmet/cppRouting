NEWS
================
Vincent LARMET
8 juillet 2019

v1.2
====

Major changes
- new functions `cpp_simplify`, `get_detour` and `to_df`

Minor changes
- bug fix in `get_distance_pair` and `get_path_pair` : verify that origin / destination nodes are present in the graph before running c++ function
- modification of `get_distance_pair` and `get_path_pair` Rd. files : clearer explanation about the choice between algorithms
- optimization of `get_distance_matrix` function : if length(to) &lt; length(from) then Dijkstra algorithm is ran from destination nodes on reversed graph
