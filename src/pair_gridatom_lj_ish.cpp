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

#include "pair_gridatom_lj_ish.hpp"

#include "neigh_list.hpp"
#include "npair.hpp"

#include <cmath>
#include <iostream>
#include <iomanip>
#include "atom.hpp"
#include "comm_brick.hpp"
#include "grid.hpp"
#include "domain.hpp"


#include "ps_pde/fftw_mpi_3darray.hpp"

using namespace PHAFD_NS;

/* ---------------------------------------------------------------------- */

PairGridAtomLJish::PairGridAtomLJish(PHAFD *phafd) : Pair(phafd) {
  list_type = NeighList::FULL;
  list_style = NPair::GRIDBIN;

  name = "gridatom/LJish";
  
};

void PairGridAtomLJish::compute()
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
  double expfac, delphi,rsq;
  
  // loop over local phi grid points which have atom neighbors

  int gridindex,jnum,aj;
  int i,j,k;

  int ajtype,ppl;

  double xexp,yexp,zexp;

  double delu,delv,delw,Ou,Ov,Ow,du,dv,dw;
  double prefac,nucwidth,ssq;

  bool printed = false;

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

	ppl = pointsperlength[ajtype];

	expfac = exp(-1./(1-rsq)+1)*epsilonstrength[ajtype];
	
	delphi = ((*grid->phi)(i,j,k)-phi[ajtype]);

	// currently no better approximation for the functional derivative of chem potential	  
	(*grid->nonlinear)(i,j,k) += 2*delphi*expfac;


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
	//       V(|\bm{r}-\bm{X}|).
	// 
	// If phi(\bm{r}) is constant, then this force should be zero by divergence theorem.
	// But this requires an exact integration of the potential V(|\bm{r}-\bm{X}|), or
	// at least a ``good enough'' numerical integration. So for us, let's opt for the latter
	// and discretise the integral by splitting space into boxes of size nucwidth**3, each
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




	if (ppl == 1) {

	  // using the value of the integrand only at the grid point


	  expfac *= 2*delphi*delphi*vol_elem/(cutsq[ajtype]*(1-rsq)*(1-rsq));
	
	  atoms->Fs(0,aj) -= delx*expfac;
	  atoms->Fs(1,aj) -= dely*expfac;
	  atoms->Fs(2,aj) -= delz*expfac;

	} else {

	  xexp = yexp = zexp = 0.0;
	  // integrating over a cube centred at the grid point, where
	  // the size of the cube is nucwidth**3
	  nucwidth = sqrt(cutsq[ajtype]);
	  Ou = delx - nucwidth/2.0;
	  Ov = dely - nucwidth/2.0;
	  Ow = delz - nucwidth/2.0;

	  du = dv = dw = nucwidth/ppl;

	  // perform integral of potential over refined grid

	  for (int iu = 0; iu < ppl; iu++) {
	    delu = iu*du + Ou;
	    for (int iv = 0; iv < ppl; iv++) {
	      delv = iv*dv + Ov;
	      for (int iw = 0; iw < ppl; iw ++) {
		delw = iw*dw + Ow;
		
		ssq = delu*delu+delv*delv+delw*delw;

		ssq /= cutsq[ajtype];


		if (ssq < 1.0) {

		  expfac = exp(-1./(1-ssq)+1)/((1-ssq)*(1-ssq));
		
		  if (commbrick->me == 3 && !printed) {
		    std::cout << "atom position = (" << atoms->xs(0,aj) 
			      << "," << atoms->xs(1,aj) << ","
			      << atoms->xs(2,aj) << ")." << std::endl;
		    std::cout << "(u,v,w) = (" << delu + atoms->xs(0,aj) << "," << delv + atoms->xs(1,aj) << "," 
			      << delw + atoms->xs(2,aj) << ")." << std::endl;
		    prefac = du*dv*dw*2*delphi*delphi/cutsq[ajtype]*epsilonstrength[ajtype];
		    std::cout << expfac*prefac*delu << std::endl;
		    std::cout << expfac*prefac*delv << std::endl;
		    std::cout << expfac*prefac*delw << std::endl;
		  }

		  xexp += expfac*delu;
		  yexp += expfac*delv;
		  zexp += expfac*delw;

		}
		
		
	      }
	    }
	  }
	  printed = true;
	  prefac = du*dv*dw*2*delphi*delphi/cutsq[ajtype]*epsilonstrength[ajtype];
	  xexp *= prefac;
	  yexp *= prefac;
	  zexp *= prefac;

	  atoms->Fs(0,aj) -= xexp;
	  atoms->Fs(1,aj) -= yexp;
	  atoms->Fs(2,aj) -= zexp;



	}
	
	// need reverse comm here to get other contributions to integrand. should be fine
	
	}

    }

  }
  return;
}


void PairGridAtomLJish::settings(const std::vector<std::string> &params)
{

  double cut = std::stod(params.at(0));
  int ppl = std::stoi(params.at(1));
  
  maxcut = cut;
  
  cutsq.resize(atoms->ntypes);
  epsilonstrength.resize(atoms->ntypes);
  phi.resize(atoms->ntypes);
  pointsperlength.resize(atoms->ntypes);


  // set some default values
  for (int i = 0; i < epsilonstrength.size(); i++) {
    cutsq[i] = cut*cut;
    epsilonstrength[i] = 1.0;
    phi[i] = 0.1;
    pointsperlength[i] = 1;
  }
  
}


void PairGridAtomLJish::coeff(const std::vector<std::string> &params)
{

  std::vector<std::string> v_line = params;
  int type = std::stoi(v_line.at(0));
  
  if (type > epsilonstrength.size() || type < 0)
    throw std::runtime_error("Pair grid/atom/LJish atom type out of range.");
  
  

  v_line.erase(v_line.begin());

  if (v_line.size() > 8) 
    throw std::runtime_error("Pair grid/atom/LJish invalid coeffs.");  
  
  int iarg = 0;
  while (iarg < v_line.size()) {

    if (v_line.at(iarg) == "epsilon") {
      epsilonstrength.at(type) = std::stod(v_line.at(iarg+1));
      iarg += 2;
    } else if (v_line.at(iarg) == "nucwidth") {
      double nucwidth = std::stod(v_line.at(iarg+1));
      cutsq.at(type) = nucwidth*nucwidth;
      if (nucwidth > maxcut) {
	throw std::runtime_error("Cannot have nucwidth greater than maxcut in grid/atom/LJish.");
      }
      iarg += 2;
    } else if (v_line.at(iarg) == "phi") {
      phi.at(type) = std::stod(v_line.at(iarg+1));
      iarg += 2;
    } else if (v_line.at(iarg) == "pointsperlength") {
      pointsperlength.at(type) = std::stoi(v_line.at(iarg+1));
      iarg += 2;
    } else
      throw std::runtime_error("Pair grid/atom/LJish invalid coeffs.");
    
  }
  
}
