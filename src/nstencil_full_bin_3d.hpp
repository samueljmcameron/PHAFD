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


#ifndef PHAFD_NSTENCIL_FULL_BIN_3D_H
#define PHAFD_NSTENCIL_FULL_BIN_3D_H

#include "nstencil.hpp"

namespace PHAFD_NS {

class NStencilFullBin3d : public NStencil {
 public:
  NStencilFullBin3d();
  void create() override;
};

}    // namespace LAMMPS_NS

#endif