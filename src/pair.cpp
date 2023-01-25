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

/* ----------------------------------------------------------------------
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "pair.hpp"
#include "neigh_list.hpp"


using namespace PHAFD_NS;

/* ---------------------------------------------------------------------- */

Pair::Pair(PHAFD *phafd) : Pointers(phafd) {
  list_type = -1;
  list_style = -1;
};


/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   specific pair style can override this function
------------------------------------------------------------------------- */

void Pair::init_list(NeighList *ptr)
{
  list = ptr;
}

