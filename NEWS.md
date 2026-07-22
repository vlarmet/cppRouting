cppRouting v4.0
===============
Major changes

-   Customizable Contraction Hierarchies (CCH) support, via `cpp_contract(..., customizable = TRUE)` and `cpp_customize()`. CCH separates preprocessing (topology) from customization (edge weights), so a prepared topology can be reused across many weight updates (e.g. congestion assignment) without a full re-contraction.
-   `get_distance_pair`, `get_distance_matrix`, `get_path_pair` and `get_aon` all support CCH graphs, in addition to normal and classically-contracted graphs.
-   `makegraph()`'s `aux` argument now accepts a `data.frame` with one or more named columns, in addition to a single vector. `aggregate_aux` in `get_distance_pair`/`get_distance_matrix` accepts `TRUE` (all columns) or a character vector of column names, and is now supported on normal, contracted (CH) and CCH graphs, and on graphs simplified with `cpp_simplify`. Aggregation runs inline in C++ during the shortest-path search itself, no path reconstruction needed.
-   `assign_traffic(..., aon_method = "cch", cch = cch)` uses a prepared CCH topology for repeated traffic assignment iterations, avoiding a full re-contraction at each run.

Minor changes

-   `cpp_simplify()` result is now classed (`cppRouting_simplified`), fixing downstream functions that silently failed on a simplified graph.
-   `cpp_simplify()` now correctly aggregates multiple auxiliary columns across merged edges (previously single-column only).
-   `cpp_contract(..., customizable = TRUE)` shows progress, like classical contraction.
-   Various CCH correctness and performance fixes: `get_path_pair` on CCH graphs, elimination ordering (near-linear instead of quadratic), many-to-many distance matrix on CCH graphs.

cppRouting v3.2
===============

- update DESCRIPTION file : remove c++11 specification and add Authors@R field

cppRouting v3.1
===============

- cppRouting use tinythread as parallel backend

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
