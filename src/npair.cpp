// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "npair.hpp"
#include <cmath>
#include "neighbor.hpp"
#include "nbin.hpp"
#include "nstencil.hpp"

#include <iostream>

#define MIN(A,B) ((A) < (B) ? (A) : (B))


using namespace PHAFD_NS;



/* ---------------------------------------------------------------------- */

NPair::NPair(PHAFD *phafd) : Pointers(phafd)
{
  last_build = -1;

}


/* ----------------------------------------------------------------------
   copy needed info from Neighbor class to this build class
   done once per run
------------------------------------------------------------------------- */

void NPair::copy_neighbor_info(const Neighbor *neighbor)
{
  // general params

  skin = neighbor->skin;
  double cutneigh = neighbor->cutneigh;
  cutneighsq = cutneigh*cutneigh;
  bboxlo = neighbor->bboxlo;
  bboxhi = neighbor->bboxhi;


}

/* ----------------------------------------------------------------------
   copy info from NBin class to this build class
------------------------------------------------------------------------- */

void NPair::copy_bin_info(const NBin *nb)
{
  nbinx = nb->nbinx;
  nbiny = nb->nbiny;
  nbinz = nb->nbinz;
  mbins = nb->mbins;
  mbinx = nb->mbinx;
  mbiny = nb->mbiny;
  mbinz = nb->mbinz;
  mbinxlo = nb->mbinxlo;
  mbinylo = nb->mbinylo;
  mbinzlo = nb->mbinzlo;

  bininvx = nb->bininvx;
  bininvy = nb->bininvy;
  bininvz = nb->bininvz;

  atom2bin = nb->atom2bin.data();
  bins = nb->bins.data();
  binhead = nb->binhead.data();
  local_atom_indices = nb->local_atom_indices;

  

}

/* ----------------------------------------------------------------------
   copy info from NStencil class to this build class
------------------------------------------------------------------------- */

void NPair::copy_stencil_info(const NStencil *ns)
{
  nstencil = ns->stencil.size();
  stencil = ns->stencil.data();

}


/* ----------------------------------------------------------------------
   copy info from NBin and NStencil classes to this build class
------------------------------------------------------------------------- */

void NPair::build_setup(const NBin *nb, const NStencil *ns,int timestep)
{
  if (nb) copy_bin_info(nb);
  if (ns) copy_stencil_info(ns);
  

  // set here, since build_setup() always called before build()
  last_build = timestep;
}

/* ----------------------------------------------------------------------
   convert grid coords into local bin #
   for orthogonal, only ghost atoms will have coord >= bboxhi or coord < bboxlo
     take special care to insure ghosts are in correct bins even w/ roundoff
     hi ghost atoms = nbin,nbin+1,etc
     owned atoms = 0 to nbin-1
     lo ghost atoms = -1,-2,etc
     this is necessary so that both procs on either side of PBC
       treat a pair of atoms straddling the PBC in a consistent way
   for triclinic, doesn't matter since stencil & neigh list built differently
------------------------------------------------------------------------- */

int NPair::coord2bin(double x, double y, double z)
{
  int ix,iy,iz;

  if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
    std::cerr << "Non-numeric positions - simulation unstable" << std::endl;


  if (x >= bboxhi[0])
    ix = static_cast<int> ((x-bboxhi[0])*bininvx) + nbinx;
  else if (x >= bboxlo[0]) {
    ix = static_cast<int> ((x-bboxlo[0])*bininvx);
    ix = MIN(ix,nbinx-1);
  } else
    ix = static_cast<int> ((x-bboxlo[0])*bininvx) - 1;

  if (y >= bboxhi[1])
    iy = static_cast<int> ((y-bboxhi[1])*bininvy) + nbiny;
  else if (y >= bboxlo[1]) {
    iy = static_cast<int> ((y-bboxlo[1])*bininvy);
    iy = MIN(iy,nbiny-1);
  } else
    iy = static_cast<int> ((y-bboxlo[1])*bininvy) - 1;

  if (z >= bboxhi[2])
    iz = static_cast<int> ((z-bboxhi[2])*bininvz) + nbinz;
  else if (z >= bboxlo[2]) {
    iz = static_cast<int> ((z-bboxlo[2])*bininvz);
    iz = MIN(iz,nbinz-1);
  } else
    iz = static_cast<int> ((z-bboxlo[2])*bininvz) - 1;

  return (iz-mbinzlo)*mbiny*mbinx + (iy-mbinylo)*mbinx + (ix-mbinxlo);
}
