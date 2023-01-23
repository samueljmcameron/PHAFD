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

#ifndef PHAFD_NPAIR_HALF_BIN_NEWTON_H
#define PHAFD_NPAIR_HALF_BIN_NEWTON_H

#include "npair.hpp"

namespace PHAFD_NS {

class NPairHalfBinNewton : public NPair {
 public:
  NPairHalfBinNewton(PHAFD *);
  void build(class NeighList *) override;

  
};

}    // namespace LAMMPS_NS

#endif
