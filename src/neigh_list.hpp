/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef PHAFD_NEIGH_LIST_H
#define PHAFD_NEIGH_LIST_H

#include "my_page.hpp"
#include <memory>
#include <vector>

namespace PHAFD_NS {

  
class NeighList {
 public:


  // data structs to store neighbor pairs I,J and associated values

  std::vector<int> ilist;          // local indices of I atoms (size nlocal)
  std::vector<int> numneigh;       // # of J neighbors for each I atom (size nlocal)
  std::vector<int *> firstneigh;    // 1st J int value of each I atom (size nlocal)

  int pgsize;            // size of each page
  int oneatom;           // max size for one atom
  std::unique_ptr<MyPage<int>> ipage;    // page of neighbor indices


  // methods

  NeighList();
  void setup_page(int, int);    // setup page data structures
};

}    // namespace LAMMPS_NS

#endif
