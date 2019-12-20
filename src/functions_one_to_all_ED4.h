#ifndef FUNCTIONS_ONE_TO_ALL_ED4_H
#define FUNCTIONS_ONE_TO_ALL_ED4_H

//Not used
void Dijkstra_mod(std::vector<std::vector<std::pair<int, double> > > &Graph,std::vector<std::vector<std::pair<int, double> > > &Graphr,
                  std::vector<std::vector<std::pair<int, double> > > &OrGraph,std::vector<std::vector<std::pair<int, double> > > &OrGraphr,
                  int dep, std::vector<int> &arr, std::vector<double> &lim, int NbNodes, std::vector<double> &Distances);

//Remove non shortest path edges in original graph
void Dijkstra_mod2(std::vector<std::vector<std::pair<int, double> > > &Graph,
                   std::vector<std::vector<std::pair<int, double> > > &OrGraph,
                   int dep, std::vector<int> &arr, std::vector<double> &lim, int NbNodes, std::vector<double> &Distances);

//Calculate number of shortcuts
int Dijkstra_bool(std::vector<std::vector<std::pair<int, double> > > &Graph,int dep, std::vector<int> &arr, std::vector<double> &lim, int NbNodes, int node, std::vector<double> &Distances);

//Edge difference
int Edge_dif(int node,std::vector<std::vector<std::pair<int, double> > > &Graph,std::vector<std::vector<std::pair<int, double> > > &Graphr,int NbNodes, std::vector<double> &Distances);

//Generate shortcuts
std::vector<std::pair<int,std::pair<int,double> > >  Dijkstra(std::vector<std::vector<std::pair<int, double> > > &Graph,
                                                              std::vector<std::vector<std::pair<int, double> > > &Graphr,
                                                              
                                                              int dep,
                                                              std::vector<int> &arr,
                                                              std::vector<double> &lim,
                                                              int NbNodes,
                                                              int node,
                                                              std::vector<double> &Distances,
                                                              std::vector<int> &Contracted,
                                                              bool reversed,
                                                              bool &err);

//Contract a node : generate shortcuts, remove contracted node
void Contract(int node,
              std::vector<std::vector<std::pair<int, double> > > &Graph,
              std::vector<std::vector<std::pair<int, double> > > &Graphr,
              std::vector<std::vector<std::pair<int, double> > > &OrGraph,
              
              int NbNodes,
              std::vector<double> &Distances,
              std::vector<int> &Contracted,
              int count,
              bool &err,
              std::vector<int> &ShortF,
              std::vector<int> &ShortT,
              std::vector<int> &ShortC);
#endif
