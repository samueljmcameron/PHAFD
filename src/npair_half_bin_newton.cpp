// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "npair_half_bin_newton.hpp"
#include "neigh_list.hpp"
#include "my_page.hpp"
#include "atom.hpp"
#include "comm_brick.hpp"

#include <iostream>



using namespace PHAFD_NS;

/* ---------------------------------------------------------------------- */

NPairHalfBinNewton::NPairHalfBinNewton(PHAFD *phafd) : NPair(phafd) {}

/* ----------------------------------------------------------------------
   binned neighbor list construction with full Newton's 3rd law
   each owned atom i checks its own bin and other bins in Newton stencil
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */


void NPairHalfBinNewton::build(NeighList *list)
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
  
  for (auto i : local_atom_indices) {

    n = 0;
    neighptr = ipage->vget();
    
    xtmp = atoms->xs(0,i);
    ytmp = atoms->xs(1,i);
    ztmp = atoms->xs(2,i);
    
    // loop over rest of atoms in i's bin, ghosts are at end of linked list
    // if j is owned atom, store it, since j is beyond i in linked list
    // if j is ghost, only store if j coords are "above and to the right" of i
    
    for (int j = bins[i-nstart]; j >= 0; j = bins[j-nstart]) {
      if (!local_atom_indices.count(j)) {
        if (atoms->xs(2,j) < ztmp) continue;
        if (atoms->xs(2,j) == ztmp) {
          if (atoms->xs(1,j) < ytmp) continue;
          if (atoms->xs(1,j) == ytmp && atoms->xs(0,j) < xtmp) continue;
        }
      }


      delx = xtmp - atoms->xs(0,j);
      dely = ytmp - atoms->xs(1,j);
      delz = ztmp - atoms->xs(2,j);
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq <= cutneighsq) {
        neighptr[n++] = j;
      }
    }

    // loop over all atoms in other bins in stencil, store every pair

    ibin = atom2bin[i-nstart];
    for (int k = 0; k < nstencil; k++) {
      for (int j = binhead[ibin+stencil[k]]; j >= 0; j = bins[j-nstart]) {
      
      
        delx = xtmp - atoms->xs(0,j);
        dely = ytmp - atoms->xs(1,j);
        delz = ztmp - atoms->xs(2,j);
        rsq = delx*delx + dely*dely + delz*delz;

        if (rsq <= cutneighsq) {
          neighptr[n++] = j;
        }
      }
    }

    ilist.push_back(i);
    firstneigh.push_back(neighptr);
    numneigh.push_back(n);
    ipage->vgot(n);
    if (ipage->status())
      std::cerr << "Neighbor list overflow, terrible things are going to happen." << std::endl;
  }


  
}
