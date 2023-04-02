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

/* ----------------------------------------------------------------------
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "pair_harmonic_cut.hpp"

#include "neigh_list.hpp"
#include "npair.hpp"

#include <cmath>


#include "atom.hpp"
#include "comm_brick.hpp"
#include "domain.hpp"


using namespace PHAFD_NS;

/* ---------------------------------------------------------------------- */

PairHarmonicCut::PairHarmonicCut(PHAFD *phafd) : Pair(phafd) {
  list_type = NeighList::HALF;
  list_style = NPair::HALFBINNEWTON;

  name = "lj/cut";
};
/* ---------------------------------------------------------------------- */



void PairHarmonicCut::compute()
{
  int i,j, jnum;
  double xtmp, ytmp, ztmp, delx, dely, delz, fpair;
  double rij;
  int *jlist;

  std::vector<int> &ilist = list->ilist;
  std::vector<int> &numneigh = list->numneigh;
  std::vector<int*> &firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (auto ii = 0; ii < ilist.size(); ii++) {

    i = ilist[ii];
    
    xtmp = atoms->xs(0,i);
    ytmp = atoms->xs(1,i);
    ztmp = atoms->xs(2,i);
    jlist = firstneigh[ii];
    jnum = numneigh[ii];



    for (int jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      delx = xtmp - atoms->xs(0,j);
      dely = ytmp - atoms->xs(1,j);
      delz = ztmp - atoms->xs(2,j);
      rij = sqrt(delx * delx + dely * dely + delz * delz);

      if (rij < cut) {

	if (rij == 0) { 

	  // if the atoms are exactly on top of each other (extremely unlikely)
	  // then just push them in a direction. should only happen when setting initial
	  // configurations...

	  delx = dely = delz = 0.001;
	  rij = sqrt(delx * delx + dely * dely + delz * delz);
	  
	}

	fpair = -1*springK*(rij-cut)/rij;


        atoms->Fs(0,i) += delx * fpair;
        atoms->Fs(1,i) += dely * fpair;
        atoms->Fs(2,i) += delz * fpair;
	
	atoms->Fs(0,j) -= delx * fpair;
	atoms->Fs(1,j) -= dely * fpair;
	atoms->Fs(2,j) -= delz * fpair;



      }
    }

  }


}

void PairHarmonicCut::settings(const std::vector<std::string>& params)
{


  cut = std::stod(params.at(0));
  springK = std::stod(params.at(1));
  maxcut = cut;



}


void PairHarmonicCut::coeff(const std::vector<std::string> & params)
{



}
