#' Compute shortest path between origin and destination nodes.
#'
#' @param Graph   An object generated by \link{makegraph}, \link{cpp_simplify} or \link{cpp_contract} function.
#' @param from A vector of one or more vertices from which shortest paths are calculated (origin).
#' @param to A vector of one or more vertices (destination).
#' @param algorithm character. \code{Dijkstra} for uni-directional Dijkstra, \code{bi} for bi-directional Dijkstra, \code{A*} for A star unidirectional search or \code{NBA} for New bi-directional A star .Default to \code{bi}
#' @param constant numeric. Constant to maintain the heuristic function admissible in A* and NBA algorithms.
#' @param keep numeric or character. Vertices of interest that will be returned.
#' @param long logical. If \code{TRUE}, a long \code{data.frame} is returned instead of a \code{list}.
#' Default to 1, when cost is expressed in the same unit than coordinates. See details
#' @return \code{list} or a \code{data.frame} containing shortest path nodes between from and to.
#' @note \code{from} and \code{from} must be the same length.
#' @details If graph is not contracted, the user has the choice between : \itemize{
#'   \item unidirectional Dijkstra (\code{Dijkstra})
#'   \item A star (\code{A*}) : projected coordinates should be provided
#'   \item bidirectional Dijkstra (\code{bi})
#'   \item New bi-directional A star (\code{NBA}) : projected coordinates should be provided
#' }
#'
#' If the input graph has been contracted by \link{cpp_contract} function, the algorithm is a modified bidirectional search.
#'
#' In \code{A*} and \code{NBA} algorithms, euclidean distance is used as heuristic function.
#'
#' All algorithms are \strong{multithreaded.} Please use \code{RcppParallel::setThreadOptions()} to set the number of threads.
#'
#' To understand the importance of constant parameter, see the package description : \url{https://github.com/vlarmet/cppRouting/blob/master/README.md}
#'
#' @seealso  \link{get_multi_paths}, \link{get_isochrone}, \link{get_detour}
#' @examples
#' #Choose number of cores used by cppRouting
#' RcppParallel::setThreadOptions(numThreads = 1)
#'
#' #Data describing edges of the graph
#' edges<-data.frame(from_vertex=c(0,0,1,1,2,2,3,4,4),
#'                   to_vertex=c(1,3,2,4,4,5,1,3,5),
#'                   cost=c(9,2,11,3,5,12,4,1,6))
#'
#' #Get all nodes
#' nodes<-unique(c(edges$from_vertex,edges$to_vertex))
#'
#' #Construct directed and undirected graph
#' directed_graph<-makegraph(edges,directed=TRUE)
#' non_directed<-makegraph(edges,directed=FALSE)
#'
#' #Sampling origin and destination nodes
#' origin<-sample(nodes,10,replace=TRUE)
#' destination<-sample(nodes,10,replace=TRUE)
#'
#' #Get distance between origin and destination in the two graphs
#' dir_paths<-get_path_pair(Graph=directed_graph, from=origin, to=destination)
#' non_dir_paths<-get_path_pair(Graph=non_directed, from=origin, to=destination)
#' print(dir_paths)
#' print(non_dir_paths)

get_path_pair<-function(Graph,from,to,algorithm="bi",constant=1,keep=NULL,long=FALSE){

  if (length(from)!=length(to)) stop("From and to have not the same length")

  if (any(is.na(cbind(from,to)))) stop("NAs are not allowed in origin/destination nodes")
  from<-as.character(from)

  to<-as.character(to)
  allnodes<-c(from,to)
  if (sum(allnodes %in% Graph$dict$ref)<length(allnodes)) stop("Some nodes are not in the graph")

  from_id<-Graph$dict$id[match(from,Graph$dict$ref)]
  to_id<-Graph$dict$id[match(to,Graph$dict$ref)]


  if (!is.null(keep)) {
    to_keep<-rep(0,Graph$nbnode)
    keep<-as.character(keep)
    to_keep[Graph$dict$ref %in% keep]<-1
  }else{
    to_keep<-rep(1,Graph$nbnode)
  }

  if (length(Graph) == 5){
    if (!is.null(Graph$coords)){
      if (algorithm %in% c("NBA","A*","bi")){
        if (algorithm=="A*"){
          if (constant == 1) warning("Are you sure constant is equal to 1 ?")
          res <- cpppath(Graph$data$from, Graph$data$to, Graph$data$dist, Graph$nbnode,
                         Graph$coords[,2],Graph$coords[,3], constant, Graph$dict$ref, to_keep, from_id, to_id, 0, 2)
        }

        if (algorithm=="NBA"){
          if (constant == 1) warning("Are you sure constant is equal to 1 ?")
          res <- cpppath(Graph$data$from, Graph$data$to, Graph$data$dist, Graph$nbnode,
                         Graph$coords[,2],Graph$coords[,3], constant, Graph$dict$ref, to_keep, from_id, to_id, 0, 3)
        }


        if (algorithm=="bi"){
          res <- cpppath(Graph$data$from, Graph$data$to, Graph$data$dist, Graph$nbnode,
                         Graph$coords[,2],Graph$coords[,3], constant, Graph$dict$ref, to_keep, from_id, to_id, 0, 1)

        }
      }

      else {
        res <- cpppath(Graph$data$from, Graph$data$to, Graph$data$dist, Graph$nbnode,
                       Graph$coords[,2],Graph$coords[,3], constant, Graph$dict$ref, to_keep, from_id, to_id, 0, 0)

      }

    }  else {
      if (algorithm=="bi"){
        res <- cpppath(Graph$data$from, Graph$data$to, Graph$data$dist, Graph$nbnode,
                       c(0,0), c(0,0), constant, Graph$dict$ref, to_keep, from_id, to_id, 0, 1)
      }
      else {
        res <- cpppath(Graph$data$from, Graph$data$to, Graph$data$dist, Graph$nbnode,
                       c(0,0), c(0,0), constant, Graph$dict$ref, to_keep, from_id, to_id, 0, 0)

      }
    }
  }


  if (length(Graph) == 6){
    res <- cpppathC(Graph$data$from, Graph$data$to, Graph$data$dist, Graph$nbnode,
                    Graph$rank, Graph$shortcuts$shortf, Graph$shortcuts$shortt, Graph$shortcuts$shortc,
                    FALSE, Graph$dict$ref, to_keep, from_id, to_id, 0)
  }



  if (long){
    names(res)<-paste0(from)
    too<-rep(to,times=sapply(res,length))
    res<-stack(setNames(res,names(res)))
    res$to<-too
    res$ind<-as.character(res$ind)
    res<-res[,c(2,3,1)]
    colnames(res)<-c("from","to","node")
    return(res)
  }else{
    names(res)<-paste0(from,"_",to)
    return(res)
  }



}
