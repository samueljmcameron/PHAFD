#include "nbin_standard.hpp"
#include "nstencil_full_bin_3d.hpp"
#include "nstencil_half_bin_3d.hpp"
#include "npair_half_bin_newton.hpp"
#include "npair_grid_bin.hpp"
#include "neighbor.hpp"
#include "neigh_list.hpp"
#include "pair.hpp"


#include "domain.hpp"
#include "atom.hpp"
#include "integrate.hpp"

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


/* must be called after pair->coeff has been created on all pairs. */
void Neighbor::setup(const std::vector<std::string> &v_line)
{

  skin = std::stod(v_line.at(0));


  double cutoff = 0.0;

  for (auto &pair : pairs) 
    cutoff = MAX(cutoff,sqrt(pair->maxcut));
  

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

  construct_lists();

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

void Neighbor::build()
{

  ago = 0;
  ncalls++;
  lastcall = integrate->timestep;

  xholds.resize(Eigen::NoChange,atoms->nowned);

  for (int i = 0; i < atoms->nowned; i++) 
    xholds.col(i) = atoms->xs.col(i); 
  
  neigh_bin->bin_atoms();

  for (auto i = 0; i < neigh_lists.size(); i ++ ) {
    neigh_pairs.at(i)->build_setup(neigh_bin.get(),neigh_stencils.at(i).get());
    neigh_pairs.at(i)->build(neigh_lists.at(i).get());
  }

  /*
  neigh_pairs[0]->build_setup(neigh_bin.get(),neigh_stencils[0].get(),timestep);
  neigh_pairs[0]->build(neigh_lists[0].get());

  neigh_pairs[1]->build_setup(neigh_bin.get(),neigh_stencils[1].get(),timestep);
  neigh_pairs[1]->build(neigh_lists[1].get());
  */
}


void Neighbor::construct_lists()
{

  
  for (auto &pair : pairs) {
    if (pair->list_type == NeighList::HALF && pair->list_style == NPair::HALFBINNEWTON) {
      bool foundlist = false;
      for (auto &nlist : neigh_lists) {
	if (nlist->list_type == NeighList::HALF && nlist->list_style == NPair::HALFBINNEWTON) {
	  pair->init_list(nlist.get());
	  foundlist = true;
	  break;
	}
      }
      if (!foundlist) {
	neigh_lists.push_back(std::make_unique<NeighList>());
	neigh_lists.back()->list_type = NeighList::HALF;
	neigh_lists.back()->list_style = NPair::HALFBINNEWTON;
	neigh_stencils.push_back(std::make_unique<NStencilHalfBin3d>());
	neigh_pairs.push_back(std::make_unique<NPairHalfBinNewton>(phafd));
	pair->init_list(neigh_lists.back().get());
      }
	
    } else if (pair->list_type == NeighList::FULL && pair->list_style == NPair::GRIDBIN) {
      bool foundlist = false;
      for (auto &nlist : neigh_lists) {
	if (nlist->list_type == NeighList::FULL && pair->list_style == NPair::GRIDBIN) {
	  pair->init_list(nlist.get());
	  foundlist = true;
	  break;
	}
      }
      if (!foundlist) {
	neigh_lists.push_back(std::make_unique<NeighList>());
	neigh_lists.back()->list_type = NeighList::FULL;
	neigh_lists.back()->list_style = NPair::GRIDBIN;
	neigh_stencils.push_back(std::make_unique<NStencilFullBin3d>());
	neigh_pairs.push_back(std::make_unique<NPairGridBin>(phafd));
	pair->init_list(neigh_lists.back().get());
      }
    } else if (pair->list_type == NeighList::NONE) {
      
      pair->init_list(nullptr);
      
    } else {
      
      throw std::runtime_error("invalid pair neighbor list requirements.");
      
    }
  }
      

  for (auto i = 0; i < neigh_lists.size(); i++) {
    neigh_lists.at(i)->setup_page(pgsize,oneatom);
    neigh_stencils.at(i)->copy_neighbor_info(this);
    neigh_stencils.at(i)->create_setup(neigh_bin.get(),integrate->timestep);
    neigh_stencils.at(i)->create();
    neigh_pairs.at(i)->copy_neighbor_info(this);
  }


  /*
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
  */

  
  return;
}
