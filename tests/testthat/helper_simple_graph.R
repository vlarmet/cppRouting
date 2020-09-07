 ## Description of Simple Graph
 ##                  +----+
 ##                  | n0 |--------+
 ##                  |    |        |9
 ##                  +----+        v
 ## +----+  12         |         +----+
 ## | n5 |<-------------------+  | n1 |
 ## |    |     +----------- --|--|    |
 ## +----+     |       |      |  +----+
 ##    |       |       |      |   ^   |
 ## 6  |       |3      |      |   |   |  11
 ##    |       |       |    +-|---+   v
 ## +-----+<---+       |    | +---+----+
 ## |n4   |<----------------|-----|n2  |
 ## |     |  5         |    |     |    |
 ## +-----+            |2   |4    +----+
 ##   |                |    |
 ##   |                v    |
 ##   |              +----+ |
 ##   + ------------>| n3 |-+
 ##         1        |    |
 ##                  +----+
 ## note shortest round n0 to n1 is n0 -> n3 -> n1 i.e. 2 + 4 = 6 i.e. shorter than direct route.
edges <- data.frame(from_vertex = c(0, 0, 1,  1, 2, 2,  3, 4, 4),
                    to_vertex =   c(1, 3, 2,  4, 4, 5,  1, 3, 5),
                    cost =        c(9, 2, 11, 3, 5, 12, 4, 1, 6))
nodes <- unique(c(edges$from_vertex, edges$to_vertex))
directed_graph <- makegraph(edges, directed = TRUE)
non_directed <- makegraph(edges, directed = FALSE)
