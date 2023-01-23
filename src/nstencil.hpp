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

#ifndef PHAFD_NSTENCIL_H
#define PHAFD_NSTENCIL_H

#include <vector>

namespace PHAFD_NS {
class NBin;
class Neighbor;

class NStencil  {
public:

  int last_stencil;    // last timestep stencil was created

  std::vector<int> stencil;                 // list of bin offsets

  int sx, sy, sz;            // extent of stencil in each dim

  NStencil();
  void copy_neighbor_info(const Neighbor *);
  virtual void create_setup(const NBin *,int);
  double memory_usage();

  virtual void create() = 0;

 protected:
  // data from Neighbor class

  double cutneigh;
  double cutneighsq;

  // data from NBin class

  int mbinx, mbiny, mbinz;
  double binsizex, binsizey, binsizez;
  double bininvx, bininvy, bininvz;


  // data common to all NStencil variants

  // methods for standard NStencil variants

  void copy_bin_info(const NBin *);                  // copy info from NBin class
  double bin_distance(int, int, int);    // distance between bin corners

};

}    // namespace LAMMPS_NS

#endif
