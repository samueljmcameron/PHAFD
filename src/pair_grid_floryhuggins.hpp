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

#ifndef PHAFD_PAIR_GRID_FLORYHUGGINS_HPP
#define PHAFD_PAIR_GRID_FLORYHUGGINS_HPP

#include "pair.hpp"

namespace psPDE {
  class FixGridFloryHuggins;
}

namespace PHAFD_NS {

  
class PairGridFloryHuggins : public Pair {
public:
  PairGridFloryHuggins(PHAFD *);
  virtual void compute() override;

  virtual void settings(const std::vector<std::string> &) override;
  virtual void coeff(const std::vector<std::string> &) override;

private:
  std::unique_ptr<psPDE::FixGridFloryHuggins> ps_flory;
  
};

}    // namespace LAMMPS_NS

#endif
