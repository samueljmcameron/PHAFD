#include "nbin_standard.hpp"
#include "nstencil_full_bin_3d.hpp"
#include "nstencil_half_bin_3d.hpp"
#include "npair_half_bin_newton.hpp"
#include "npair_grid_bin.hpp"
#include "neighbor.hpp"
#include "neigh_list.hpp"

#include "domain.hpp"
#include "atom.hpp"

#include <iostream>

#define MAX(A,B) ((A) > (B) ? (A) : (B))

using namespace PHAFD_NS;

Neighbor::Neighbor(PHAFD *phafd) : Pointers(phafd)
{
  binsizeflag = false;
  binsize_user = 0.0;
  cutneigh = 0.0;
  pgsize = 100000;
  oneatom = 2000;
  skin = 0.0;

  delay = 0;
  every = 1;
  dist_check = true;
  ndanger = 0;
  ncalls = 0;
}

Neighbor::~Neighbor() = default;

void Neighbor::setup(double cutoff, double skintmp)
{
  skin = skintmp;
  triggersq = 0.25*skin*skin;
  
  if (cutoff > 0.0) {
    cutneigh = skin+cutoff;
  }

  
  
  // build neighbor bins, and also copy in box dimensions
  neigh_bin = std::make_unique<NBinStandard>(phafd);
  bboxlo = domain->boxlo;
  bboxhi = domain->boxhi;
  neigh_bin->copy_neighbor_info(this);
  neigh_bin->setup_bins(cutoff);

  // construct a stencil for a full neighbor and half neighbor list,
  //  store full stencil in [0], half stencil in [1]
  neigh_stencils.push_back(std::make_unique<NStencilFullBin3d>());
  neigh_stencils.back()->copy_neighbor_info(this);
  neigh_stencils.back()->create_setup(neigh_bin.get(),0);
  neigh_stencils.back()->create();
  neigh_stencils.push_back(std::make_unique<NStencilHalfBin3d>());
  neigh_stencils.back()->copy_neighbor_info(this);
  neigh_stencils.back()->create_setup(neigh_bin.get(),0);
  neigh_stencils.back()->create();


  // two lists, one for full and one for half
  neigh_lists.push_back(std::make_unique<NeighList>());
  neigh_lists.back()->setup_page(pgsize,oneatom);
  neigh_lists.push_back(std::make_unique<NeighList>());
  neigh_lists.back()->setup_page(pgsize,oneatom);

  // two neighbor pair setups, one for full, one for half
  neigh_pairs.push_back(std::make_unique<NPairGridBin>(phafd));
  neigh_pairs.back()->copy_neighbor_info(this);  
  neigh_pairs.push_back(std::make_unique<NPairHalfBinNewton>(phafd));
  neigh_pairs.back()->copy_neighbor_info(this);


}

bool Neighbor::decide()
{
  ago++;
  if (ago >= delay && ago % every == 0) {
    if (!dist_check) return true;
    return check_distance();
  } return false;
}

bool Neighbor::check_distance()
{

  int flag = 0;
  double rsq;
  for (int i = 0; i < atoms->nowned; i++) {
    rsq = (atoms->xs.col(i)-xholds.col(i)).squaredNorm();
    if (rsq > triggersq) flag = 1;
  }

  int flagall;

  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_MAX,world);
  if (flagall && ago == MAX(every,delay)) ndanger++;
  if (flagall) return true;
  return false;
  

}

void Neighbor::build(int timestep)
{

  ago = 0;
  ncalls++;
  lastcall = timestep;

  xholds.resize(Eigen::NoChange,atoms->nowned);

  for (int i = 0; i < atoms->nowned; i++) 
    xholds.col(i) = atoms->xs.col(i); 
  
  neigh_bin->bin_atoms(timestep);

  neigh_pairs[0]->build_setup(neigh_bin.get(),neigh_stencils[0].get(),timestep);
  neigh_pairs[0]->build(neigh_lists[0].get());

  neigh_pairs[1]->build_setup(neigh_bin.get(),neigh_stencils[1].get(),timestep);
  neigh_pairs[1]->build(neigh_lists[1].get());
}

