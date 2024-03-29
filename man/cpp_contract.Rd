% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/contract.R
\name{cpp_contract}
\alias{cpp_contract}
\title{Contraction hierarchies algorithm}
\usage{
cpp_contract(Graph, silent = FALSE)
}
\arguments{
\item{Graph}{An object generated by \link{makegraph} or \link{cpp_simplify} function.}

\item{silent}{Logical. If \code{TRUE}, progress is not displayed.}
}
\value{
A contracted graph.
}
\description{
Contract a graph by using contraction hierarchies algorithm
}
\details{
Contraction hierarchies is a speed-up technique for finding shortest path in a graph.
It consist of two steps : preprocessing phase and query. \code{cpp_contract()} preprocess the input graph to later use special query algorithm implemented in \link{get_distance_pair}, \link{get_distance_matrix}, \link{get_aon} and \link{get_path_pair} functions.
To see the benefits of using contraction hierarchies, see the package description : \url{https://github.com/vlarmet/cppRouting/blob/master/README.md}.
}
\examples{
#Data describing edges of the graph
edges<-data.frame(from_vertex=c(0,0,1,1,2,2,3,4,4),
                  to_vertex=c(1,3,2,4,4,5,1,3,5),
                  cost=c(9,2,11,3,5,12,4,1,6))

#Construct cppRouting graph
graph<-makegraph(edges,directed=TRUE)

#Contract graph
contracted_graph<-cpp_contract(graph,silent=TRUE)
}
\seealso{
\link{cpp_simplify}
}
