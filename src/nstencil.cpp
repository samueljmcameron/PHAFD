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

#include "nstencil.hpp"
#include "neighbor.hpp"
#include "nbin.hpp"

using namespace PHAFD_NS;

/* ----------------------------------------------------------------------
   NStencil classes
   each has method to create a stencil = list of bin offsets
     invoked each time simulation box size/shape changes
     since induces change in bins
   stencil = bins whose closest corner to central bin is within cutoff
   sx,sy,sz = bin bounds = furthest the stencil could possibly extend
     calculated below in create_setup()
   3d creates xyz stencil, 2d creates xy stencil
   for full list or half list with newton off
     use a full stencil
     stencil is all surrounding bins including self
     regardless of triclinic
   for half list with newton on
     use a half stencil
     stencil is bins to the "upper right" of central bin
     stencil does not include self
     no versions that allow ghost on (no callers need it?)
   for half list with newton on and triclinic:
     use a half stencil
     stencil is all bins in z-plane of self and above, but not below
     in 2d is all bins in y-plane of self and above, but not below
     stencil includes self
     no versions that allow ghost on (no callers need it?)
   for multi/old:
     create one stencil for each atom type
     stencil follows same rules for half/full, newton on/off, triclinic
     cutoff is not cutneighmaxsq, but max cutoff for that atom type
     no versions that allow ghost on (any need for it?)
   for multi:
     create one stencil for each icollection-jcollection pairing
     the full/half stencil label refers to the same-collection stencil
     a half list with newton on has a half same-collection stencil
     a full list or half list with newton off has a full same-collection stencil
     cross collection stencils are always full to allow small-to-large lookups
     for orthogonal boxes, a half stencil includes bins to the "upper right" of central bin
     for triclinic, a half stencil includes bins in the z (3D) or y (2D) plane of self and above
     cutoff is not cutneighmaxsq, but max cutoff for that atom collection
     no versions that allow ghost on (any need for it?)
------------------------------------------------------------------------- */

NStencil::NStencil()
{
  last_stencil = -1;

}


/* ----------------------------------------------------------------------
   copy needed info from Neighbor class to this stencil class
------------------------------------------------------------------------- */

void NStencil::copy_neighbor_info(const Neighbor * neighbor)
{
  cutneigh = neighbor->cutneigh;
  cutneighsq = cutneigh*cutneigh;
}

/* ----------------------------------------------------------------------
   copy needed info from NBin class to this stencil class
------------------------------------------------------------------------- */

void NStencil::copy_bin_info(const NBin * nb)
{
  mbinx = nb->mbinx;
  mbiny = nb->mbiny;
  mbinz = nb->mbinz;
  binsizex = nb->binsizex;
  binsizey = nb->binsizey;
  binsizez = nb->binsizez;
  bininvx = nb->bininvx;
  bininvy = nb->bininvy;
  bininvz = nb->bininvz;
}

/* ----------------------------------------------------------------------
   insure NBin data is current
   insure stencils are allocated large enough
------------------------------------------------------------------------- */

void NStencil::create_setup(const NBin *nb,int timestep)
{


  if (nb) copy_bin_info(nb);
  last_stencil = timestep;

  // sx,sy,sz = max range of stencil in each dim
  // smax = max possible size of entire 3d stencil
  // stencil will be empty if cutneigh = 0.0
  
  sx = static_cast<int>(cutneigh*bininvx);
  if (sx*binsizex < cutneigh) sx++;
  sy = static_cast<int>(cutneigh*bininvy);
  if (sy*binsizey < cutneigh) sy++;
  sz = static_cast<int>(cutneigh*bininvz);
  if (sz*binsizez < cutneigh) sz++;

  
}

/* ----------------------------------------------------------------------
   compute closest distance between central bin (0,0,0) and bin (i,j,k)
------------------------------------------------------------------------- */

double NStencil::bin_distance(int i, int j, int k)
{
  double delx,dely,delz;

  if (i > 0) delx = (i-1)*binsizex;
  else if (i == 0) delx = 0.0;
  else delx = (i+1)*binsizex;

  if (j > 0) dely = (j-1)*binsizey;
  else if (j == 0) dely = 0.0;
  else dely = (j+1)*binsizey;

  if (k > 0) delz = (k-1)*binsizez;
  else if (k == 0) delz = 0.0;
  else delz = (k+1)*binsizez;

  return (delx*delx + dely*dely + delz*delz);
}
