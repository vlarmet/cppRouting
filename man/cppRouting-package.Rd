\name{cppRouting-package}
\alias{cppRouting}
\docType{package}
\title{
  Fast Implementation of Dijkstra Algorithm in R
}
\description{
  Calculation of distances, shortest paths and isochrones on non-negative weighted graphs using several variants of Dijkstra algorithm.
}
\section{Functions}{
  
\itemize{
\item distance matrix (between all combinations origin-destination nodes),
\item distances between origin and destination by pair (one-to-one),
\item shortest paths between origin and destination by pair (one-to-one),
\item shortest paths between all origin nodes and all destination nodes,
\item Isochrones/isodistances with one or multiple breaks.
\item graph simplification by removing non intersection nodes, duplicated edges and isolated loops
\item find nodes that can be reached under an additional detour time around the shortest path between two nodes
}
}

\section{Algorithms}{
  
  Algorithms can be chosen for one-to-one calculations like get_distance_pair() and get_path_pair() :
\itemize{
\item uni-directional Dijkstra algorithm,
\item bi-directional Dijkstra algorithm,
\item uni-directional A* algorithm (nodes coordinates are needed),
\item New bi-directional A* algorithm (nodes coordinates are needed).
}
}


\references{
Dijkstra, E. W. (1959), A note on two problems in connexion with graphs.

P. E. Hart, N. J. Nilsson et B. Raphael (1968). A Formal Basis for the Heuristic Determination of Minimum Cost Paths.

Goldberg, Andrew & Fonseca F. Werneck, Renato (2005). Computing Point-to-Point Shortest Paths from External Memory.

Pijls & Post (2009). Yet another bidirectional algorithm for shortest paths.

}



