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

#include "domain.hpp"
#include "atom.hpp"
#include "ps_pde/grid.hpp"
#include <vector>
#include <string>

namespace PHAFD {
  
class NeighList;
  
class Pair {
public:
  
  NeighList *list;        // standard neighbor list used by most pairs

  Pair(Atom *, psPDE::Grid *);




  // general child-class methods

  virtual void compute(const Domain &) = 0;
  
  virtual void init_list(NeighList *ptr);

  virtual void coeff(std::vector<std::string> ) = 0;

  virtual void settings(std::vector<std::string>) = 0;

protected:
  Atom *atoms;

  psPDE::Grid *grid;

};

}    // namespace LAMMPS_NS

#endif
