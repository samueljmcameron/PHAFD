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

/* ----------------------------------------------------------------------
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "pair_gridatom_lj_ish_linear.hpp"

#include "neigh_list.hpp"
#include "npair.hpp"

#include <cmath>
#include <iostream>
#include <iomanip>
#include "atom.hpp"
#include "comm_brick.hpp"
#include "grid.hpp"
#include "domain.hpp"


#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;

/* ---------------------------------------------------------------------- */

PairGridAtomLJishLinear::PairGridAtomLJishLinear(PHAFD *phafd) : Pair(phafd) {
  list_type = NeighList::FULL;
  list_style = NPair::GRIDBIN;

  name = "gridatom/LJish/linear";

};

void PairGridAtomLJishLinear::compute()
{

  double xtmp,ytmp,ztmp,delx,dely,delz;
  int *jlist;

  std::vector<int> &ilist = list->ilist;
  std::vector<int> &numneigh = list->numneigh;
  std::vector<int*> &firstneigh = list->firstneigh;

  const double Ox = domain->sublo[0];
  const double Oy = domain->sublo[1];
  const double Oz = domain->sublo[2];

  const double vol_elem = grid->dx()*grid->dy()*grid->dz();
  double expfac,rsq;
  
  // loop over local phi grid points which have atom neighbors

  int gridindex,jnum,aj;
  int i,j,k;

  int ajtype;




  std::cout << std::fixed << std::setprecision(12);  
  for (int iibin = 0; iibin < ilist.size(); iibin++) { 

    gridindex = ilist[iibin];

    grid->phi->reverseFlat(gridindex,i,j,k);

    xtmp = k*grid->dx()+Ox;
    ytmp = j*grid->dy()+Oy;
    ztmp = i*grid->dz()+Oz;

    jlist = firstneigh[iibin];
    jnum = numneigh[iibin];

    // then over atom neighbors
    for (int ajj = 0; ajj < jnum; ajj++) { 
      aj = jlist[ajj];

      ajtype = atoms->types[aj];
    
      delx = xtmp - atoms->xs(0,aj);
      dely = ytmp - atoms->xs(1,aj);
      delz = ztmp - atoms->xs(2,aj);

      rsq = delx*delx + dely*dely + delz*delz;
	    
      if (rsq < cutsq[ajtype]) {

	// rescale rsq by the nucleation width of the atom
	rsq /= cutsq[ajtype];

	expfac = -1*exp(-1./(1-rsq)+1)*epsilonstrength[ajtype];

	// currently no better approximation for the functional derivative of chem potential	  
	(*grid->chempot)(i,j,k) += expfac;


	// better approximation for forces felt by atoms interacting with grid is necessary
	// so that lattice effects are minimal
       
	// current thinking:
	// free energy integral is
	//
	//      U = \int d^3 r (\phi(\bm{r})-\bar{\phi})^2 V(|\bm{r}-\bm{X}|),
	//
	// where \bm{X} is the position of the atom (single atom here for simpilcity).
	// The (local) chemical potential is
	//
	//  \frac{\delta U}{\delta \phi} = 2(\phi(\bm{r})-\bar{\phi})V(|\bm{r}-\bm{X}|),
	//
	// and the force on the atom is
	//
	// -\frac{\partial U}{\partial \bm{X}} 
	//  = - \int d^3 r (\phi(\bm{r})-\bar{\phi})^2 \frac{\partial}{\partial \bm{X}}
	//       V(|\bm{r}-\bm{X}|)
	//
	// or by changing the integration variable to u = \bm{r} -\bm{X} and back again,
	// -\frac{\partial U}{\partial \bm{X}} 
	//  = - \int d^3 r 2(\phi(\bm{r})-\bar{\phi})\nabla\phi(\bm{r}) V(|\bm{r}-\bm{X}|)
	// 
	// If phi(\bm{r}) is constant, then this force should be zero either by divergence theorem
	// (for first version of force above), or zero explicitly by second version above.
	// First version requires an exact integration of the potential V(|\bm{r}-\bm{X}|), or
	// at least a ``good enough'' numerical integration so the divergence theorem approximately
	// holds (otherwise, the particle will experience significant forces when \phi\neq \bar{\phi}
	// entirely due to numerical error). So for us, let's opt for the latter version of the force,
	// which requires an extra round of FFTW'ing to compute gradphi, but often we need to do this
	// anyway so it should be fine.
	// we discretise the integral by splitting space into boxes of size nucwidth**3, each
	// box being centred on a grid point. Each box is then further discretised into 
	// pointsperlength**3 points, where the potential is evaluated at these points but
	// the grid is held constant at phi(\bm{r}), where \bm{r} is the grid point that the
	// box is centred around. This isn't perfect by any means, but for large enough
	// pointsperlength, the result that the force is zero for constant phi will be recovered
	// at least. So the resulting integral has a contribution from each box of the form
	//
	// 
	//  -(\phi(\bm{r})-\bar{\phi})^2 du*dv*dw sum_{u,v,w}(  \frac{\partial}{\partial \bm{X}}
	//       V(|\bm{s}-\bm{X}|).
	//
	// where u\in[-nucwidth/2.0+x,nucwidth/2.0+x) (and similar for v and w) and
	// \bm{s} = (u,v,w) are the points to be evaluated in the integral, and 
	// \bm{r} = (x,y,z) is the grid point (at the centre of the u,v,w domain).



	
	// using the value of the integrand only at the grid point
	
	
	expfac *= vol_elem;
	
	atoms->Fs(0,aj) -= (*grid->gradphi[0])(i,j,k)*expfac;
	atoms->Fs(1,aj) -= (*grid->gradphi[1])(i,j,k)*expfac;
	atoms->Fs(2,aj) -= (*grid->gradphi[2])(i,j,k)*expfac;
	
	
	// need reverse comm here to get other contributions to integrand. should be fine
	
      }
      
    }
    
  }

  return;
}


void PairGridAtomLJishLinear::settings(const std::vector<std::string> &params)
{

  double cut = std::stod(params.at(0));
  
  maxcut = cut;
  
  cutsq.resize(atoms->ntypes);
  epsilonstrength.resize(atoms->ntypes);



  // set some default values
  for (int i = 0; i < epsilonstrength.size(); i++) {
    cutsq[i] = cut*cut;
    epsilonstrength[i] = 1.0;
  }
  
}


void PairGridAtomLJishLinear::coeff(const std::vector<std::string> &params)
{

  std::vector<std::string> v_line = params;
  int type = std::stoi(v_line.at(0));
  
  if (type > epsilonstrength.size() || type < 0)
    throw std::runtime_error("Pair grid/atom/LJish/linear atom type out of range.");
  
  

  v_line.erase(v_line.begin());

  if (v_line.size() > 6) 
    throw std::runtime_error("Pair grid/atom/LJish/linear invalid coeffs.");  
  
  int iarg = 0;
  while (iarg < v_line.size()) {

    if (v_line.at(iarg) == "epsilon") {
      epsilonstrength.at(type) = std::stod(v_line.at(iarg+1));
      iarg += 2;
    } else if (v_line.at(iarg) == "nucwidth") {
      double nucwidth = std::stod(v_line.at(iarg+1));
      cutsq.at(type) = nucwidth*nucwidth;
      if (nucwidth > maxcut) {
	throw std::runtime_error("Cannot have nucwidth greater than maxcut in grid/atom/LJish/linear.");
      }
      iarg += 2;
    } else
      throw std::runtime_error("Pair grid/atom/LJish/linear invalid coeffs.");
    
  }
  
}
