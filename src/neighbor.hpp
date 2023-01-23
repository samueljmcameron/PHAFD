

#ifndef PHAFD_NEIGHBOR_HPP
#define PHAFD_NEIGHBOR_HPP


#include <Eigen/Core>

#include <memory>
#include <array>

#include "pointers.hpp"

namespace PHAFD_NS {
class NBin;
class NStencil;
class NeighList;
class NPair;
  
class Neighbor : protected Pointers {
  
  int delay,every,ago,ncalls,lastcall;
  bool dist_check;
  double triggersq;

  bool check_distance();
  
  
public:
  Neighbor(PHAFD *);
  ~Neighbor();
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
  const double *bboxlo,*bboxhi;
  void setup(double,double);
  void build(int);
  bool decide();

  int get_ncalls() {return ncalls;};

};
  
}
#endif
