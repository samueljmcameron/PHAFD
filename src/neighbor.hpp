

#ifndef PHAFD_NEIGHBOR_HPP
#define PHAFD_NEIGHBOR_HPP


#include <Eigen/Core>

#include "domain.hpp"
#include "atom.hpp"
#include "ps_pde/grid.hpp"

#include <memory>
#include <array>
#include <mpi.h>



namespace PHAFD {
class NBin;
class NStencil;
class NeighList;
class NPair;
  
class Neighbor {
  
  int delay,every,ago,ncalls,lastcall;
  bool dist_check;
  double triggersq;

  bool check_distance(Atom &, MPI_Comm);
  
  
public:
  Neighbor();
  Eigen::Matrix3Xd xholds;

  int ndanger;
  bool binsizeflag;
  double binsize_user;
  double cutneigh,skin;
  int pgsize;
  int oneatom;
  std::unique_ptr<NBin> neigh_bin;
  std::vector<std::unique_ptr<NStencil>> neigh_stencils;
  std::vector<std::unique_ptr<NeighList>> neigh_lists;
  std::vector<std::unique_ptr<NPair>> neigh_pairs;
  std::array<double,3> bboxlo,bboxhi;
  void setup(const psPDE::Domain &, double,double);
  void build(const Atom &,const psPDE::Grid &,int);
  bool decide(Atom &, MPI_Comm);

  int get_ncalls() {return ncalls;};

};
  
}
#endif
