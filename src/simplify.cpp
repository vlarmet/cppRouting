#include "graph.h"
#include <vector>
#include <Rcpp.h>

using namespace std;

void quickDelete( int idx ,vector<pair<int,double> > &vec)
{
  vec[idx] = vec.back();
  vec.pop_back();
}

void quickDelete2( int idx ,vector<int> &vec)
{
  vec[idx] = vec.back();
  vec.pop_back();
}

void Graph::simp(vector<IVec> &Gr,
                 IVec &Junction,
                 vector<IVec > &Edges,
                 IVec &keep,
                 IVec &Treated,
                 bool &loop){
  // Detect junction nodes
  for (int i=0; i < data.size() ; ++i){
    //Nodes to keep
    if (keep[i]==1){
      continue;
    }

    //undirected
    if ((data[i].size()==2 && Gr[i].size()==2)){
      vector<pair<int,double> > vec1=data[i];
      IVec vec2 = Gr[i];
      sort(vec1.begin(),vec1.end());
      sort(vec2.begin(),vec2.end());

      if (vec1[0].first==vec2[0] && vec1[1].first==vec2[1]){
        Junction[i]=2;

      }


    }

    //directed
    if (data[i].size()==1 && Gr[i].size()==1){
      if (data[i][0].first!=Gr[i][0]){
        Junction[i]=1;


      }


    }

  }
  ////////////////////////////////////////////////////////////////////////////////
  for (int i=0; i < data.size() ; ++i){
    if (Junction[i]!=2 || Treated[i]==1){
      continue;
    }



    //Rcpp::Rcerr << i << endl;
    Rcpp::checkUserInterrupt();

    IVec path;
    int node1=i;
    Treated[i]=1;
    path.push_back(i);




    if (Gr[node1].size()!=2){
      continue;
    }

    int breaked=0;

    for (int j=0; j<2; ++j){
      int node2=Gr[node1][j];

      breaked=0;
      path.insert(path.begin(),node2);
      Treated[node2]=1;

      while (Junction[node2]==2){
        Rcpp::checkUserInterrupt();


        int node=Gr[node2][0];


        if (Treated[node]==1 && Junction[node]==2){

          node=Gr[node2][1];
          //loop in the data
          if (Treated[node]==1 && Junction[node]==2){
            breaked=1;
            break;
          }

        }
        //Rcpp::Rcerr << node << endl;

        node2=node;

        //loop in the data

        if (node==i){
          breaked=1;
          break;
        }

        path.insert(path.begin(),node2);
        Treated[node2]=1;

      }

      reverse(path.begin(),path.end());


    }



    // for (int k=0; k < path.size();++k){
    //  Rcpp::Rcerr << path[k] << "__";
    //}


    if (breaked==0){


      Edges.push_back(path);
    }




  }


  int breaked=0;
  /// Construct directed edges
  for (int i=0; i < data.size() ; ++i){
    if (Junction[i]!=1 || Treated[i]==1){
      continue;
    }

    breaked=0;


    //Rcpp::Rcerr << i << endl;
    Rcpp::checkUserInterrupt();

    IVec path;
    int node1=i;
    path.push_back(i);
    Treated[i]=1;

    if (Gr[node1].size()!=1){
      continue;
    }

    int node2=Gr[node1][0];

    path.insert(path.begin(),node2);
    Treated[node2]=1;

    //One direction
    while (Junction[node2]==1){
      Rcpp::checkUserInterrupt();

      int node=Gr[node2][0];

      //Rcpp::Rcerr << node << endl;

      node2=node;
      //loop in the data
      if (node==i){
        breaked=1;
        break;
      }


      path.insert(path.begin(),node2);
      Treated[node2]=1;

    }

    //other direction

    node2=data[node1][0].first;
    path.push_back(node2);
    Treated[node2]=1;
    while (Junction[node2]==1){
      Rcpp::checkUserInterrupt();

      int node=data[node2][0].first;

      //Rcpp::Rcerr << node << endl;

      node2=node;
      //loop in the data
      if (node==i){
        breaked=1;
        break;
      }


      path.push_back(node2);
      Treated[node2]=1;

    }


    if (breaked==0){
      Edges.push_back(path);
    }




  }

  /////////////////////////////////////////////////////////////////////////////////
  IVec Removenode(data.size(),0);  // nodes inserted in new edge


  DVec Acc(Edges.size());
  DVec Acc2(Edges.size());
  for (int i=0; i < Edges.size(); ++i){
    IVec path=Edges[i];
    double w=0.0;
    for (int k=0; k < (path.size()-1); ++k){

      Removenode[path[k]]=1;
      int node1=path[k];
      int node2=path[k+1];

      for (int m=0; m < data[node1].size(); ++m){
        if (data[node1][m].first==node2){
          w+=data[node1][m].second;
        }
      }

    }




    Acc[i]=w;
    //Autre sens
    IVec path2=path;
    reverse(path2.begin(),path2.end());
    double w2=0.0;
    for (int k=0; k < (path2.size()-1); ++k){
      int node1=path2[k];
      int node2=path2[k+1];

      for (int m=0; m < data[node1].size(); ++m){
        if (data[node1][m].first==node2){
          w2+=data[node1][m].second;
        }
      }



    }
    Acc2[i]=w2;

  }

  ///////////////////////////////////////////////////////////////////////
  //Remove from data and Gr
  for (int i=0; i < Edges.size();i++){
    for (int j=1; j < Edges[i].size()-1;j++){
      //data
      vector<pair<int,double > > Gr_remov=data[Edges[i][j]];

      data[Edges[i][j]].erase(data[Edges[i][j]].begin(),data[Edges[i][j]].end());
      for (int k=0; k < Gr[Edges[i][j]].size();k++){
        for (int m=0; m < data[Gr[Edges[i][j]][k]].size();m++){
          if (data[Gr[Edges[i][j]][k]][m].first==Edges[i][j]) {
            quickDelete(m,data[Gr[Edges[i][j]][k]]);
          }
        }

      }

      //Gr
      Gr[Edges[i][j]].erase(Gr[Edges[i][j]].begin(),Gr[Edges[i][j]].end());
      for (int k=0; k < Gr_remov.size();k++){
        for (int m=0; m < Gr[Gr_remov[k].first].size();m++){
          if (Gr[Gr_remov[k].first][m]==Edges[i][j]){

            quickDelete2(m,Gr[Gr_remov[k].first]);
          }
        }
      }
    }
  }


  //Remove loops
  if (loop==true){
    for (int i=0; i < data.size(); ++i){
      if (Junction[i]>0 && Removenode[i]==0){
        vector<pair<int,double > > Gr_remov=data[i];

        data[i].erase(data[i].begin(),data[i].end());
        for (int k=0; k < Gr[i].size();k++){
          for (int m=0; m < data[Gr[i][k]].size();m++){
            if (data[Gr[i][k]][m].first==i) {
              quickDelete(m,data[Gr[i][k]]);
            }
          }

        }

        Gr[i].erase(Gr[i].begin(),Gr[i].end());
        for (int k=0; k < Gr_remov.size();k++){
          for (int m=0; m < Gr[Gr_remov[k].first].size();m++){
            if (Gr[Gr_remov[k].first][m]==i){

              quickDelete2(m,Gr[Gr_remov[k].first]);
            }
          }
        }
      }

    }

  }

  //Add shortcuts and avoid duplicated edges
  bool dup;
  for (int i=0; i < Edges.size(); i++){
    if (Edges[i][0]==Edges[i][Edges[i].size()-1] && keep[Edges[i][0]]==0){
      continue;
    }


    dup=false;
    int node1=Edges[i][0];
    int node2=Edges[i][Edges[i].size()-1];

    for (int k=0; k<data[node1].size(); k++){
      if (data[node1][k].first==node2){
        dup=true;
        if (data[node1][k].second >= Acc[i]){
          data[node1][k].second=Acc[i];
          break;
        }


      }
    }

    if (dup==false) {
      data[node1].push_back(make_pair(node2,Acc[i]));
      Gr[node2].push_back(node1);
    }

    //Other direction
    if (Junction[Edges[i][1]]==2){
      dup=false;
      for (int k=0; k<data[node2].size(); k++){
        if (data[node2][k].first==node1){
          dup=true;
          if (data[node2][k].second >= Acc2[i]){
            data[node2][k].second=Acc2[i];
            break;
          }


        }
      }

      if (dup==false) {
        data[node2].push_back(make_pair(node1,Acc2[i]));
        Gr[node1].push_back(node2);
      }
    }


  }


  ///////////////////////////////////////////////////////////////////
  //Number of nodes
  IVec Nodes(data.size(),0);
  for (int i=0; i < data.size();i++){
    for (int j=0; j < data[i].size();j++){
      Nodes[data[i][j].first]+=1;
    }

    for (int j=0; j < Gr[i].size();j++){
      Nodes[Gr[i][j]]+=1;
    }

  }

  int count=0;
  for (int i=0; i < Nodes.size();i++){
    if (Nodes[i]>0) count+=1;
  }

  nbnode = count;

  ////////////////////////////////////////////////////////////////
  //Reinitialize
  fill(Treated.begin(),Treated.end(),0);
  fill(Junction.begin(),Junction.end(),0);
  Edges.resize(0);
}



///////////////////////////////////////////////////////////////////

void Graph::simplify(bool loop,IVec keep,
                     bool iterate, bool progress){



  std::vector<IVec> Gr(nbnode);


  for (unsigned int i = 0; i < nbnode; ++i) {

    for (unsigned int j = 0; j < data[i].size(); j++){
      Gr[data[i][j].first].push_back(i);
    }

  }

  IVec Treated(nbnode,0);
  IVec Junction(nbnode,0);
  std::vector<IVec> Edges;


  int count=0;

  if (iterate){
    while (true){

      int N2 = nbnode;
      simp(Gr, Junction, Edges, keep, Treated, loop);
      count+=1;
      if (progress && (N2 - nbnode)!=0) Rcpp::Rcout << "iteration : "<<count<< " - " <<N2-nbnode<< " nodes removed" << std::endl;
      if (N2 == nbnode) break;

    }

  }
  else{
    int N2 = nbnode;
    simp(Gr,Junction,Edges,keep,Treated,loop);

    if (progress) Rcpp::Rcout<< N2-nbnode << " nodes removed" << std::endl;
  }


  // Update number of edges
  count=0;
  for (int i = 0; i < nbnode; i++) count += data[i].size();
  nbedge = count;

}
