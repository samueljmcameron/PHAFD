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

#ifndef PHAFD_PAIR_LJ_CUT_H
#define PHAFD_PAIR_LJ_CUT_H

#include "pair.hpp"

namespace PHAFD {

class PairLJCut : public Pair {
public:
  PairLJCut(Atom *, psPDE::Grid *);
  void compute(const Domain &) override;

  void settings( std::vector<std::string> ) override;
  void coeff( std::vector<std::string> ) override;

private:
  double lj1,lj2;

  double cutsq;
  
  
};

}    // namespace LAMMPS_NS

#endif
