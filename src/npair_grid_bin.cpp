

#include "npair_grid_bin.hpp"
#include "neigh_list.hpp"
#include "my_page.hpp"

#include "atom.hpp"
#include "grid.hpp"

#include "ps_pde/fftw_mpi_3darray.hpp"

#include <iostream>

using namespace PHAFD_NS;

NPairGridBin::NPairGridBin(PHAFD *phafd) : NPair(phafd) {}

/* ----------------------------------------------------------------------
   binned neighbor list construction for all neighbors of every grid point
------------------------------------------------------------------------- */

void NPairGridBin::build(NeighList *list)
{
  int n,ibin;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr;

  
  
  
  std::vector<int> &ilist = list->ilist;
  std::vector<int> &numneigh = list->numneigh;

  
  std::vector<int*> &firstneigh = list->firstneigh;
  
  MyPage<int> *ipage = list->ipage.get();

  ilist.clear();
  numneigh.clear();
  firstneigh.clear();
  ipage->reset();

  int nstart = atoms->nowned;
  
  int totalgridbins = grid->phi->Nx()*grid->phi->Ny()*grid->phi->Nz();
  

  /*
   MyPage<int> is a data structure which allocates large spaces at a time for storing
    index values of neighboring atoms, e.g. when grid point 0 has 3 atom neighbors,
    grid point 1 has 2 atom neighbors, grid point 2 has 0 neighbors, grid point 3
    has 2 neighbors, grid point 4 has 6 neighbors, and MyPage has a page size of 10,
    then the following values will be written for the relevant output:
    inum = 5
    ilist = {0,1,2,3,4}
    firstneigh = {pointer to page0, index 0 ; pointer to page0, index 3 ;
                  pointer to page0, index 5 ; pointer to page0, index 5,
		  pointer to page1, index 0 }

    ipage page 1 = {g0s first neighbor, g0s second neighbor, g0s thirdneighbor,
                    g1s first neighbor, g2s first neighbor, g3s first neighbor,
	            g3s second neighbor, notallocd, notalloc, notallocd}
    ipage page 2 = {g4s first neighbor, g4s second neighbor, ... g4s sixth neighbor,
                    notallocd,notallocd,notallocd,notallocd}

    numneigh = {3,2,0,2,6}

  */

  int gridindex = 0;

  int local0start = grid->phi->get_local0start();
  for (int k = 0; k < grid->phi->Nz(); k++) {

    ztmp = (k+local0start)*grid->dz()+bboxlo[2];

    for (int j = 0; j < grid->phi->Ny(); j++) {

      ytmp = j*grid->dy()+bboxlo[1];
      
      for (int i = 0; i < grid->phi->Nx(); i++) {

	n = 0;
	neighptr = ipage->vget();


	xtmp = i*grid->dx()+bboxlo[0];

	
	ibin = coord2bin(xtmp,ytmp,ztmp);
	

	for (int m = 0; m < nstencil; m++) {
	  for (int ai = binhead[ibin + stencil[m]]; ai >= 0; ai = bins[ai-nstart]) {
	    
	    delx = xtmp - atoms->xs(0,ai);
	    dely = ytmp - atoms->xs(1,ai);
	    delz = ztmp - atoms->xs(2,ai);
	    
	    rsq = delx*delx + dely*dely + delz*delz;
	    
	    if (rsq <= cutneighsq) {
	      neighptr[n++] = ai;
	    }
	  }
	}

	if (n > 0) {
	
	  ilist.push_back(gridindex++);
	  firstneigh.push_back(neighptr);
	  numneigh.push_back(n);
	} else {
	  gridindex ++;
	}
	ipage->vgot(n);
	
	if (ipage->status())
	  std::cout << "Neighbor list overflow, boost neigh_modify one" << std::endl;

      }
    }
  }

}
