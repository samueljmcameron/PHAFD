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


#ifndef PHAFD_NBIN_STANDARD_H
#define PHAFD_NBIN_STANDARD_H

#include "nbin.hpp"

namespace PHAFD {

class NBinStandard : public NBin {
 public:
  NBinStandard();

  void setup_bins(const psPDE::Domain &, double) override;
  void bin_atoms(const Atom &, int) override;

};

}    // namespace LAMMPS_NS

#endif
