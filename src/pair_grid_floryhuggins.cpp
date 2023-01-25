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

#include "pair_grid_floryhuggins.hpp"

#include "neigh_list.hpp"

#include <cmath>

#include "atom.hpp"
#include "comm_brick.hpp"
#include "grid.hpp"
#include "domain.hpp"

#include "ps_pde/fixgrid_floryhuggins.hpp"

using namespace PHAFD_NS;

/* ---------------------------------------------------------------------- */

PairGridFloryHuggins::PairGridFloryHuggins(PHAFD *phafd) : Pair(phafd) {

  ps_flory = std::make_unique<psPDE::FixGridFloryHuggins>();
  
};

void PairGridFloryHuggins::compute()
{


  ps_flory->compute(*(grid->ps_grid));
  
}


void PairGridFloryHuggins::settings(const std::vector<std::string> &params)
{

  ps_flory->readCoeffs(params);

}


void PairGridFloryHuggins::coeff(const std::vector<std::string> &params)
{

}
