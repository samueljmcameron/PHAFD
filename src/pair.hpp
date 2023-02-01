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

#ifndef PHAFD_PAIR_H
#define PHAFD_PAIR_H

#include <vector>
#include <string>

#include "pointers.hpp"

namespace PHAFD_NS {
  
class NeighList;
  
class Pair : protected Pointers {
public:
  
  NeighList *list;        // standard neighbor list used by most pairs

  Pair(PHAFD *);

  // general child-class methods

  virtual void compute() = 0;
  
  virtual void init_list(NeighList *ptr);

  virtual void coeff(const std::vector<std::string> &) = 0;

  virtual void settings(const std::vector<std::string> &) = 0;

  double maxcut;

  int list_type,list_style;
  std::string name;

protected:

};

}    // namespace LAMMPS_NS

#endif
