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

#ifndef PHAFD_NPAIR_H
#define PHAFD_NPAIR_H

#include <set>
#include "atom.hpp"
#include "ps_pde/grid.hpp"

namespace PHAFD {
  
class Neighbor;
class NBin;
class NStencil;
class NeighList;


  
class NPair {
 public:


  int last_build;     // last timestep build performed

  NPair();
  virtual void copy_neighbor_info(const Neighbor *);
  void build_setup(const NBin *, const NStencil *, int);
  virtual void build(class NeighList *,const  Atom &) = 0;
  virtual void build(class NeighList *, const Atom &,
		     const psPDE::Grid &) = 0;

 protected:

  // data from Neighbor class

  double skin;
  double cutneighsq;
  const double *bboxlo,*bboxhi;


  // data from NBin class

  int nbinx, nbiny, nbinz;
  int mbins;
  int mbinx, mbiny, mbinz;
  int mbinxlo, mbinylo, mbinzlo;
  double bininvx, bininvy, bininvz;
  const int *binhead;     // index of first atom in each bin
  const int *bins;        // index of next atom in same bin
  const int *atom2bin;    // bin assignment for each atom (local+ghost)
  std::set<int> local_atom_indices; // indices of local atoms
  

  // data from NStencil class

  int nstencil;
  const int *stencil;


  // methods for all NPair variants

  virtual void copy_bin_info(const NBin *);
  virtual void copy_stencil_info(const NStencil *);

  int coord2bin(double,double,double);


};

}    // namespace LAMMPS_NS

#endif
