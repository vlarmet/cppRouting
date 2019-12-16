#ifndef SIMP_H
#define SIMP_H

void Simplify(std::vector<std::vector<std::pair<int,double> > > &Graph,
              std::vector<std::vector<int> > &Gr,
              std::vector<int> &Junction,
              std::vector<std::vector<int> > &Edges,
              std::vector<int> &keep,
              std::vector<int> &Treated,
              bool &loop,
              int &N);
#endif