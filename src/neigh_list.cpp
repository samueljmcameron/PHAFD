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

#include "neigh_list.hpp"
#include "my_page.hpp"


using namespace PHAFD_NS;

#define PGDELTA 1

/* ---------------------------------------------------------------------- */

NeighList::NeighList() : ipage(std::make_unique<MyPage<int>>())
{


}

/* ---------------------------------------------------------------------- */

void NeighList::setup_page(int pgsize_caller, int oneatom_caller)
{
  pgsize = pgsize_caller;
  oneatom = oneatom_caller;

  ipage->init(oneatom,pgsize,PGDELTA);


}

