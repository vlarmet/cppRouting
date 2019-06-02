cppRouting package
================
Vincent LARMET
30 mai 2019

Package presentation
====================

`cppRouting` is a `R` package which provide functions to calculate distances, shortest paths and isochrones/isodistances on weighted graphs. Under the hood, `cppRouting` call Dijkstra and A\* algorithms implemented in C++ using std::priority\_queue from the Standard Template Library.
This package have been made with `Rcpp` and `RcppParallel`.

Install from github
===================

``` r
library(devtools)
devtools::install_github("vlarmet/cppRouting")
```

Data
====

The data presented here is the official french road network describing over 500000 km of roads.
All data used in this README are free and can be downloaded here :

-   roads : <http://professionnels.ign.fr/route500>
-   general practitioners location : <https://www.insee.fr/fr/statistiques/3568614?sommaire=3568656#consulter>
-   maternity wards location : <https://www.insee.fr/fr/statistiques/3568611?sommaire=3568656#dictionnaire>
-   shapefile of the ~36000 communes in France : <http://professionnels.ign.fr/adminexpress>

Graph data have been preprocessed for more readability (see data\_preparation.R).

The final graph is composed of 234615 nodes and 685118 edges.
Data has to be a 3 columns data.frame or matrix containing from, to and a cost/distance column. Here the cost is the time needed to travel in each edges (in minutes). From and to are vertices IDs (character or numeric).

Main functions
==============

`cppRouting` package provide these functions :
- distance matrix (between all combinations origin-destination nodes),
- distances between origin and destination by pair,
- shortest paths between origin and destination by pair,
- shortest paths between all origin nodes and all destination nodes,
- Isochrones/isodistances with one or multiple breaks.

The choice between Dijkstra and A\* algorithm is available for `get_distance_pair` and `get_path_pair`. In these functions, Dijkstra algorithm is stopped when the destination node is reached.
A\* is relevant if geographic coordinates of all nodes are provided. Note that coordinates should be expressed in a projection system.
To be accurate and efficient, `A*` algorithm should use an admissible heuristic function (see <https://en.wikipedia.org/wiki/A*_search_algorithm>), e.g the cost and geographic coordinates must be expressed in the same unit.
In `cppRouting`, heuristic function `h` is defined such that : h(xi,yi,xdestination,ydestination)/k, with a constant k; so in the case where coordinates are expressed in meters and cost is expressed in time, k is the maximum speed allowed on the road.

By default, constant is 1 and is designed for graphs with cost expressed in the same unit than coordinates (e.g meters).

### Let's see the benefit of the A\* algorithm with the french road network :

``` r
library(cppRouting)
library(dplyr)
library(sf)
library(ggplot2)
#setwd("")
#Reading french road data
roads<-read.csv("roads.csv",colClasses = c("character","character","numeric"))
#Shapefile data of communes
com<-read_sf("com_simplified_geom.shp")
#Correspondance file between communes and nodes in the graph
ndcom<-read.csv("node_commune.csv",colClasses = c("character","character","numeric"))
#General practitioners locations
med<-read.csv("doctor.csv",colClasses = c("character","numeric","character","numeric"))
#Import nodes coordinates (projected in EPSG : 2154)
coord<-read.csv("coordinates.csv",colClasses = c("character","numeric","numeric"))
#Head of road network data
head(roads)
```

    ##   from     to    weight
    ## 1    0 224073 0.4028571
    ## 2    1  65036 3.5280000
    ## 3    2 173723 1.8480000
    ## 4    3      2 2.5440000
    ## 5    4 113129 4.9680000
    ## 6    5      4 1.6680000

### Head of coordinates data

``` r
head(coord)
```

    ##   ID        X       Y
    ## 1  0 805442.8 6458384
    ## 2  1 552065.9 6790520
    ## 3  2 556840.2 6790475
    ## 4  3 554883.7 6790020
    ## 5  4 548345.2 6791000
    ## 6  5 547141.3 6790434

### Instantiate the graph

``` r
#Instantiate a graph with coordinates
graph<-makegraph(roads,directed = T,coords = coord)
```

### Run Dijkstra algorithm for finding minimum cost between pairs of nodes

``` r
#Generate 2000 random origin and destination nodes
origin<-sample(roads$from,2000)
destination<-sample(roads$from,2000)

#Benchmarks single core
#Dijkstra 
system.time(
pair_dijkstra<-get_distance_pair(graph,origin,destination)
)
```

    ## Running Dijkstra ...

    ##    user  system elapsed 
    ##   54.69    0.00   54.79

``` r
#Benchmarks parallel
#Dijkstra
system.time(
pair_dijkstra_par<-get_distance_pair(graph,origin,destination,allcores = TRUE)
)
```

    ## Running Dijkstra ...

    ##    user  system elapsed 
    ##   70.31    0.00   17.97

### Run A\* algorithm

Coordinates are defined in meters and max speed is 110km/h; so for the heuristic function to be admissible, the constant equal 110/0.06 :

``` r
#A* single node
system.time(
pair_astar<-get_distance_pair(graph,origin,destination,algorithm = "A*",constant = 110/0.06)
)
```

    ## Running A* ...

    ##    user  system elapsed 
    ##   30.42    2.28   32.79

``` r
#A* parallel
system.time(
pair_astar_par<-get_distance_pair(graph,origin,destination,algorithm = "A*",constant = 110/0.06,allcores = TRUE)
)
```

    ## Running A* ...

    ##    user  system elapsed 
    ##   44.38    0.72   11.65

A\* is the fastest one and the output is the same.

``` r
head(cbind(pair_dijkstra,pair_astar,pair_dijkstra_par,pair_astar_par))
```

    ##      pair_dijkstra pair_astar pair_dijkstra_par pair_astar_par
    ## [1,]      462.4514   462.4514          462.4514       462.4514
    ## [2,]      114.1351   114.1351          114.1351       114.1351
    ## [3,]      238.4923   238.4923          238.4923       238.4923
    ## [4,]      490.4173   490.4173          490.4173       490.4173
    ## [5,]      218.5715   218.5715          218.5715       218.5715
    ## [6,]      482.1238   482.1238          482.1238       482.1238

Applications
============

Application 1 : Calculate Two Step Floating Catchment Areas (2SFCA) of general practitioners in France
------------------------------------------------------------------------------------------------------

2SFCA method is explained here : <https://en.wikipedia.org/wiki/Two-step_floating_catchment_area_method>

### First step

Isochrones are calculated with the `cppRouting` function `get_isochrone`

``` r
#Isochrone around doctor locations with time limit of 15 minutes
iso<-get_isochrone(graph,from = ndcom[ndcom$com %in% med$CODGEO,"id_noeud"],lim = 15)
#Convert list to long data frame
df<-stack(setNames(iso, seq_along(iso)))
df$ind<-rep(names(iso),times=sapply(iso,length))
df<-df[df$values %in% ndcom$id_noeud,]
#Joining and summing population located in each isochrone
df<-left_join(df,ndcom[,c("id_noeud","POPULATION")],by=c("values"="id_noeud"))
df<-df %>% group_by(ind) %>%
  summarise(pop=sum(POPULATION))
#Joining number of doctors 
df<-left_join(df,med[,c("id_noeud","NB_D201")],by=c("ind"="id_noeud"))
#Calculate ratios
df$ratio<-df$NB_D201/df$pop
```

### Second step

``` r
#Isochrone around each commune with time limit of 15 minutes (few seconds to compute)
iso2<-get_isochrone(graph,from=ndcom$id_noeud,lim = 15)
#Convert list to long data frame
df2<-stack(setNames(iso2, seq_along(iso2)))
df2$ind<-rep(names(iso2),times=sapply(iso2,length))
#Joining and summing ratios calculated in first step
df2<-left_join(df2,df[,c("ind","ratio")],by=c("values"="ind"))
df2<-df2 %>% group_by(ind) %>%
  summarise(sfca=sum(ratio,na.rm=T))
```

### Plot the map for Bourgogne-Franche-Comte region

``` r
#Joining commune IDs to nodes
df2<-left_join(df2,ndcom[,c("id_noeud","com")],by=c("ind"="id_noeud"))
#Joining 2SFCA to shapefile
com<-left_join(com,df2[,c("com","sfca")],by=c("INSEE_COM"="com"))
#Plot for one region
p<-ggplot()+
  geom_sf(data=com[com$NOM_REG=="BOURGOGNE-FRANCHE-COMTE",],aes(fill=sfca),colour=NA)+
  coord_sf(datum=NA)+
  scale_fill_gradient(low="#BB2528",high = "#FFFF66")+
  labs(fill="2SFCA")+
  ggtitle("2SFCA applied to general practitioners")
p
```

![](readme_files/figure-markdown_github/unnamed-chunk-8-1.png)

Application 2 : Calculate the minimum travel time to the closest maternity ward in France
-----------------------------------------------------------------------------------------

``` r
#Import materinty ward locations
maternity<-read.csv("maternity.csv",colClasses = c("character","numeric"))
```

### Shortest travel time matrix

The shortest travel time is computed with the `cppRouting` function `get_distance_matrix`.
We compute travel time from all commune nodes to all maternity ward nodes (e.g ~36000\*400 distances).

``` r
#Distance matrix (around 10 minutes to compute)
dists<-get_distance_matrix(graph,
                           from=ndcom$id_noeud,
                           to=ndcom$id_noeud[ndcom$com %in% maternity$CODGEO],
                           allcores=TRUE)
#We extract each minimum travel time for all the communes
dists2<-data.frame(node=ndcom$id_noeud,mindist=apply(dists,1,min,na.rm=T))
#Joining commune IDs to nodes
dists2<-left_join(dists2,ndcom[,c("id_noeud","com")],by=c("node"="id_noeud"))
#Joining minimum travel time to the shapefile
com<-left_join(com,dists2[,c("com","mindist")],by=c("INSEE_COM"="com"))
```

Plot the map of minimum travel time in Bourgogne-Franche-Comte region
=====================================================================

``` r
p<-ggplot()+
  geom_sf(data=com[com$NOM_REG=="BOURGOGNE-FRANCHE-COMTE",],aes(fill=mindist),colour=NA)+
  coord_sf(datum=NA)+
  scale_fill_gradient(low="#009900",high="#003300")+
  labs(fill="Minutes")+
  ggtitle("Travel time to the closest maternity ward")
p
```

![](readme_files/figure-markdown_github/unnamed-chunk-11-1.png)

Benchmark with other R packages
===============================

To show the efficiency of `cppRouting`, we can make some benchmarking with the famous R package `igraph`, and the `dodgr` package which provide highly optimized heaps.

### Distance matrix : one core

``` r
library(igraph)
library(dodgr)
#Sampling 1000 random origin/destination nodes (1000000 distances to compute)
origin<-sample(unique(roads$from),1000,replace = F)
destination<-sample(unique(roads$from),1000,replace = F)
```

``` r
#igraph 
graph_igraph<-graph_from_data_frame(roads,directed = TRUE)

system.time(
  test_igraph<-distances(graph_igraph,origin,to=destination,weights = E(graph_igraph)$weight,mode="out")
)
```

    ##    user  system elapsed 
    ##   86.48    0.06   86.69

``` r
#dodgr
#Adding coordinates to data
roads2<-roads
colnames(roads2)[3]<-"dist"
roads2<-left_join(roads2,coord,by=c("from"="ID"))
colnames(roads2)[4:5]<-c("from_lon","from_lat")
roads2<-left_join(roads2,coord,by=c("to"="ID"))
colnames(roads2)[6:7]<-c("to_lon","to_lat")
colnames(roads2)[1:2]<-c("from_id","to_id")
roads2$from_id<-as.character(roads2$from_id)
roads2$to_id<-as.character(roads2$to_id)

system.time(
test_dodgr<-dodgr_dists(graph=data.frame(roads2),from=origin,to=destination,parallel=FALSE)
)
```

    ##    user  system elapsed 
    ##   85.02    0.08   85.27

``` r
#cppRouting
system.time(
test_cpp<-get_distance_matrix(graph,origin,destination,allcores = FALSE)
)
```

    ##    user  system elapsed 
    ##   54.35    0.02   54.45

### Distance matrix : parallel

``` r
#dodgr
system.time(
test_dodgr<-dodgr_dists(graph=data.frame(roads2),from=origin,to=destination,parallel=TRUE)
)
```

    ##    user  system elapsed 
    ##  118.82    0.36   31.42

``` r
#cppRouting
system.time(
test_cpp<-get_distance_matrix(graph,origin,destination,allcores = TRUE)
)
```

    ##    user  system elapsed 
    ##   69.43    0.00   17.73

Benchmarking on shortest paths by pairs
---------------------------------------

``` r
#Sampling 500 random origin/destination nodes 
origin<-sample(unique(roads$from),500,replace = F)
destination<-sample(unique(roads$from),500,replace = F)
#dodgr
system.time(
test_dodgr<-dodgr_paths(graph=data.frame(roads2),from=origin,to=destination,pairwise = TRUE)
)
```

    ##    user  system elapsed 
    ##  524.96   18.16  544.13

``` r
#cppRouting
system.time(
test_cpp<-get_path_pair(graph,origin,destination,algorithm = "A*",constant=110/0.06)
)
```

    ## Running A* ...

    ##    user  system elapsed 
    ##    7.88    0.01    7.90

### Test similarity of the first travel

``` r
#Number of nodes
length(test_dodgr[[1]][[1]])
```

    ## [1] 324

``` r
length(test_cpp[[1]])
```

    ## [1] 324

``` r
#Setdiff 
setdiff(test_dodgr[[1]][[1]],test_cpp[[1]])
```

    ## character(0)

New functions `cppRouting` will provide in the future
=====================================================

-   Detours admitting shortest paths : finding the nodes that are reachable under a fixed detour time around the shortest path
-   Graph simplification by removing irrelevant nodes in order to compute in a faster way the shortest distance or travel time
-   Contraction hierarchies implementation
