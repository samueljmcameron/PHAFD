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

#include "pair_lj_cut.hpp"

#include "neigh_list.hpp"
#include "npair.hpp"

#include <cmath>


#include "atom.hpp"
#include "comm_brick.hpp"
#include "domain.hpp"



using namespace PHAFD_NS;

/* ---------------------------------------------------------------------- */

PairLJCut::PairLJCut(PHAFD *phafd) : Pair(phafd) {
  list_type = NeighList::HALF;
  list_style = NPair::HALFBINNEWTON;

  name = "lj/cut";
};
/* ---------------------------------------------------------------------- */



void PairLJCut::compute()
{
  int i,j, jnum;
  double xtmp, ytmp, ztmp, delx, dely, delz, fpair;
  double rsq, r2inv, r6inv, forcelj, factor_lj;
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
      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cutsq) {
        r2inv = 1.0 / rsq;
        r6inv = r2inv * r2inv * r2inv;
        forcelj = r6inv * (lj1*r6inv - lj2);
        fpair = forcelj * r2inv;

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

void PairLJCut::settings(const std::vector<std::string>& params)
{


  double cut = std::stod(params.at(0));
  double epsilon = std::stod(params.at(1));
  double sigma = std::stod(params.at(2));

  maxcut = cut;
  
  lj1 = 48*epsilon*pow(sigma,12.0);
  lj2 = 24*epsilon*pow(sigma,6.0);
  cutsq = cut*cut;



}


void PairLJCut::coeff(const std::vector<std::string> & params)
{



}
