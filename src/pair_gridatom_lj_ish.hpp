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

#ifndef PHAFD_PAIR_GRIDATOM_LJISH_HPP
#define PHAFD_PAIR_GRIDATOM_LJISH_HPP

#include "pair.hpp"

namespace PHAFD_NS {

class PairGridAtomLJish : public Pair {
public:
  PairGridAtomLJish(PHAFD *);
  virtual void compute() override;

  virtual void settings(const std::vector<std::string> &) override;
  virtual void coeff(const std::vector<std::string> &) override;

private:
  std::vector<double> epsilonstrength,phi,cutsq;

};

}    // namespace LAMMPS_NS

#endif
