/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
p   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef PHAFD_NBIN_H
#define PHAFD_NBIN_H

#include <set>
#include <vector>

#include "pointers.hpp"

namespace PHAFD_NS {

class Neighbor;
  
class NBin : protected Pointers {
 public:
  int last_bin;         // last timestep atoms were binned

  // Variables for NBinStandard

  int nbinx, nbiny, nbinz;    // # of global bins
  int mbins;                  // # of local bins and offset on this proc
  int mbinx, mbiny, mbinz;
  int mbinxlo, mbinylo, mbinzlo;

  std::set<int> local_atom_indices; // indices of local atoms
  // (to be used for iterating over in npair).

  double binsizex, binsizey, binsizez;    // bin sizes and inverse sizes
  double bininvx, bininvy, bininvz;

  std::vector<int> binhead;     // index of first atom in each bin (size is total # of bins)
  std::vector<int> bins;        // index of next atom in same bin (size is ngathered)
  std::vector<int> atom2bin;    // bin assignment for each atom (local+ghost) (size is ngathered)

  NBin(PHAFD *);
  virtual void copy_neighbor_info(const Neighbor *);

  virtual void setup_bins(double) = 0;
  virtual void bin_atoms(int) = 0;

  

 protected:
  // data from Neighbor class

  double cutneigh;
  int binsizeflag;
  double binsize_user;
  const double *bboxlo, *bboxhi;

  // methods

  int coord2bin(const double *);

};

}    // namespace LAMMPS_NS

#endif
