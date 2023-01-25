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

#include <algorithm>

#include "atom.hpp"
#include "nbin_standard.hpp"
#include "domain.hpp"
#include "neighbor.hpp"

using namespace PHAFD_NS;

#define SMALL 1.0e-6
#define CUT2BIN_RATIO 100

/* ---------------------------------------------------------------------- */

NBinStandard::NBinStandard(PHAFD *phafd) : NBin(phafd) {}

/* ----------------------------------------------------------------------
   setup neighbor binning geometry
   bin numbering in each dimension is global:
     0 = 0.0 to binsize, 1 = binsize to 2*binsize, etc
     nbin-1,nbin,etc = bbox-binsize to bbox, bbox to bbox+binsize, etc
     -1,-2,etc = -binsize to 0.0, -2*binsize to -binsize, etc
   code will work for any binsize
     since next(xyz) and stencil extend as far as necessary
     binsize = 1/2 of cutoff is roughly optimal
   for orthogonal boxes:
     a dim must be filled exactly by integer # of bins
     in periodic, procs on both sides of PBC must see same bin boundary
     in non-periodic, coord2bin() still assumes this by use of nbin xyz
   for triclinic boxes:
     tilted simulation box cannot contain integer # of bins
     stencil & neigh list built differently to account for this
   mbinlo = lowest global bin any of my ghost atoms could fall into
   mbinhi = highest global bin any of my ghost atoms could fall into
   mbin = number of bins I need in a dimension
------------------------------------------------------------------------- */

void NBinStandard::setup_bins(double cutoff)
{
  // bbox = size of bbox of entire domain
  // bsubbox lo/hi = bounding box of my subdomain extended by comm->cutghost
  // for triclinic:
  //   bbox bounds all 8 corners of tilted box
  //   subdomain is in lamda coords
  //   include dimension-dependent extension via comm->cutghost
  //   domain->bbox() converts lamda extent to box coords and computes bbox

  std::array<double,3> bbox,bsubboxlo,bsubboxhi;
  double cutghost = cutoff;


  bsubboxlo = {domain->sublo[0],domain->sublo[1],domain->sublo[2]};
  bsubboxlo[0] -= cutghost;
  bsubboxlo[1] -= cutghost;
  bsubboxlo[2] -= cutghost;

  bsubboxhi = {domain->subhi[0],domain->subhi[1],domain->subhi[2]};
  bsubboxhi[0] += cutghost;
  bsubboxhi[1] += cutghost;
  bsubboxhi[2] += cutghost;


  bbox[0] = bboxhi[0] - bboxlo[0];
  bbox[1] = bboxhi[1] - bboxlo[1];
  bbox[2] = bboxhi[2] - bboxlo[2];

  // optimal bin size is roughly 1/2 the cutoff
  // for BIN style, binsize = 1/2 of max neighbor cutoff
  // for MULTI_OLD style, binsize = 1/2 of min neighbor cutoff
  // special case of all cutoffs = 0.0, binsize = box size

  double binsize_optimal;
  if (binsizeflag) binsize_optimal = binsize_user;
  else binsize_optimal = 0.5*cutneigh;
  if (binsize_optimal == 0.0) binsize_optimal = bbox[0];
  
  double binsizeinv = 1.0/binsize_optimal;


  // create actual bins
  // always have one bin even if cutoff > bbox
  // for 2d, nbinz = 1

  nbinx = static_cast<int> (bbox[0]*binsizeinv);
  nbiny = static_cast<int> (bbox[1]*binsizeinv);
  nbinz = static_cast<int> (bbox[2]*binsizeinv);


  if (nbinx == 0) nbinx = 1;
  if (nbiny == 0) nbiny = 1;
  if (nbinz == 0) nbinz = 1;

  // compute actual bin size for nbins to fit into box exactly
  // error if actual bin size << cutoff, since will create a zillion bins
  // this happens when nbin = 1 and box size << cutoff
  // typically due to non-periodic, flat system in a particular dim
  // in that extreme case, should use NSQ not BIN neighbor style

  binsizex = bbox[0]/nbinx;
  binsizey = bbox[1]/nbiny;
  binsizez = bbox[2]/nbinz;

  bininvx = 1.0 / binsizex;
  bininvy = 1.0 / binsizey;
  bininvz = 1.0 / binsizez;

  if (binsize_optimal*bininvx > CUT2BIN_RATIO ||
      binsize_optimal*bininvy > CUT2BIN_RATIO ||
      binsize_optimal*bininvz > CUT2BIN_RATIO)
    throw std::runtime_error("Cannot use neighbor bins - box size << cutoff");

  // mbinlo/hi = lowest and highest global bins my ghost atoms could be in
  // coord = lowest and highest values of coords for my ghost atoms
  // static_cast(-1.5) = -1, so subract additional -1
  // add in SMALL for round-off safety

  int mbinxhi,mbinyhi,mbinzhi;
  double coord;

  coord = bsubboxlo[0] - SMALL*bbox[0];
  mbinxlo = static_cast<int> ((coord-bboxlo[0])*bininvx);
  if (coord < bboxlo[0]) mbinxlo = mbinxlo - 1;
  coord = bsubboxhi[0] + SMALL*bbox[0];
  mbinxhi = static_cast<int> ((coord-bboxlo[0])*bininvx);

  coord = bsubboxlo[1] - SMALL*bbox[1];
  mbinylo = static_cast<int> ((coord-bboxlo[1])*bininvy);
  if (coord < bboxlo[1]) mbinylo = mbinylo - 1;
  coord = bsubboxhi[1] + SMALL*bbox[1];
  mbinyhi = static_cast<int> ((coord-bboxlo[1])*bininvy);


  coord = bsubboxlo[2] - SMALL*bbox[2];
  mbinzlo = static_cast<int> ((coord-bboxlo[2])*bininvz);
  if (coord < bboxlo[2]) mbinzlo = mbinzlo - 1;
  coord = bsubboxhi[2] + SMALL*bbox[2];
  mbinzhi = static_cast<int> ((coord-bboxlo[2])*bininvz);
  

  // extend bins by 1 to insure stencil extent is included
  // for 2d, only 1 bin in z

  mbinxlo = mbinxlo - 1;
  mbinxhi = mbinxhi + 1;
  mbinx = mbinxhi - mbinxlo + 1;

  mbinylo = mbinylo - 1;
  mbinyhi = mbinyhi + 1;
  mbiny = mbinyhi - mbinylo + 1;


  mbinzlo = mbinzlo - 1;
  mbinzhi = mbinzhi + 1;
  mbinz = mbinzhi - mbinzlo + 1;

  int bbin = ((int) mbinx) * ((int) mbiny) * ((int) mbinz) + 1;
  mbins = bbin;

  binhead.resize(mbins);
}

/* ----------------------------------------------------------------------
   bin owned and ghost atoms
------------------------------------------------------------------------- */

void NBinStandard::bin_atoms(int step)
{
  int ibin;


  
  last_bin = step;

  for (auto &b: binhead) b = -1;

  local_atom_indices.clear();


  // bin in reverse order so linked list will be in forward order
  // also puts ghost atoms at end of list, which is necessary

  int nstart = atoms->nowned;
  int nall = atoms->nowned + atoms->ngathered;

  atom2bin.resize(nall-nstart);
  bins.resize(nall-nstart);

  // loop through ghost atoms first, storing local atom indices in an (ordered) set
  //  to use for later
  for (int i = nall-1; i >= nstart; i--) {
    if (atoms->labels[i] == Atom::GHOST) {
      ibin = coord2bin(&(atoms->xs(0,i)));
      atom2bin[i-nstart] = ibin;
      bins[i-nstart] = binhead[ibin];
      binhead[ibin] = i;
    } else { // local_atom_indices is a set, so will order indices from small to large
      local_atom_indices.insert(i);
    }
  }

  // now loop over local atom indices, which are in order from smallest to largest
  
  std::set<int>::reverse_iterator rit;

  for (rit=local_atom_indices.rbegin(); rit != local_atom_indices.rend(); ++rit) {
    ibin = coord2bin(&(atoms->xs(0,*rit)));
    atom2bin[*rit-nstart] = ibin;
    bins[*rit-nstart] = binhead[ibin];
    binhead[ibin] = *rit;
  }

  
  
}

