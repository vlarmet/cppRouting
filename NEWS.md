cppRouting v3.0
===============
Major changes

-   All C++ code have been rewritten in a much more OOP way, and can be easily used outside of Rcpp
-   All routing algorithms are now implemented in subclasses of RcppParallel::Worker and are natively multithreaded
-   graph are now internally represented as adjacency lists (3 vectors), which is generally 20% faster for routing
-   `get_path_pair`, `get_isochrone`, `get_detour` and `get_multi_paths` are now multithreaded
-   aggregating a secondary weight along shortest path is now implemented for `get_distance_pair` and `get_distance_matrix`, for both normal and contracted network, and multithreaded
-   implementation of all-or-nothing assignment for both normal and contracted network, multithreaded
-   link-based algorithms for the calculation of the User Equilibrium (UE) : Method of Successive Averages, variants of Frank-Wolfe algorithm (normal, conjugate, bi-conjugate). Multithreaded
-   bush-based algorithm for the calculation of the User Equilibrium (UE) : algorithm-B from R.B. Dial (with batching for large OD matrix). Partially multithreaded.
-   implementation of stall-on-demand technique (see Geisberger, 2008) which speed-up routing computation on a contracted network (both pairwise and matrix)

Minor changes

-   bug fix in `get_path_pair` and `get_multi_paths` : node sequence listed from `origin` to `destination` node.


cppRouting v2.0
===============

Major changes

-   implementation of contraction hierarchies algorithm
-   thread-safe implementation of all parallel algorithms with `RcppParallel` package instead of `parallel`
-   implementation of one-to-one query algorithm on contracted graph 
-   implementation of many-to-many query algorithm on contracted graph 
-   implementation of PHAST algorithm on contracted graph 

Minor changes

-   new options `long` and `keep` for `get_path_pair`, `get_isochrone` and `get_multi_paths` functions    
-   remove `allcores` option for `get_detour` function  
-   optimization of `cpp_simplify` function  
-   optimization of `makegraph` function

cppRouting v1.2
===============

Major changes

-   new functions `cpp_simplify`, `get_detour` and `to_df`

Minor changes

-   bug fix in `get_distance_pair` and `get_path_pair` : verify that origin / destination nodes are present in the graph before running c++ function
-   modification of `get_distance_pair` and `get_path_pair` Rd. files : clearer explanation about the choice between algorithms
-   optimization of `get_distance_matrix` function : if length(to) &lt; length(from) then Dijkstra algorithm is ran from destination nodes on reversed graph

cppRouting v1.1
===============

First CRAN release
