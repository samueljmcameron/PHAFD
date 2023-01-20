#include "pair_gridatom_gaussian.hpp"

#include "neigh_list.hpp"
#include <cmath>
#include <iostream>

using namespace PHAFD;

/* ---------------------------------------------------------------------- */

PairGridAtomGaussian::PairGridAtomGaussian(Atom *atoms,psPDE::Grid *grid)
  : Pair(atoms,grid)
{}



// interaction between field and polymer
// goes in a pair_*.cpp file

void PairGridAtomGaussian::compute(const Domain & domain)
{



  double xtmp,ytmp,ztmp,delx,dely,delz;
  int *jlist;

  std::vector<int> &ilist = list->ilist;
  std::vector<int> &numneigh = list->numneigh;
  std::vector<int*> &firstneigh = list->firstneigh;

  const double Ox = domain.sublo[0];
  const double Oy = domain.sublo[1];
  const double Oz = domain.sublo[2];

  const double dv = grid->dx()*grid->dy()*grid->dz();
  double expfac, delphi,rsq;
  
  // loop over local phi grid points which have atom neighbors

  int gridindex,jnum,aj;
  int i,j,k;

  int ajtype;
  
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
	    
      if (rsq <= cutsq[ajtype]) {

	// rescale rsq by the nucleation width of the atom
	rsq /= nucwidth[ajtype]*nucwidth[ajtype];
	
	expfac = exp(-rsq/2.0)*epsilonstrength[ajtype];
	
	delphi = ((*grid->phi)(i,j,k)-phi[ajtype]);

	(*grid->nonlinear)(i,j,k) += 2*delphi*expfac;

	// re-use expfac vs creating new variable
	expfac *= delphi*delphi*dv;
	
	atoms->Fs(0,aj) -= delx*expfac;
	atoms->Fs(1,aj) -= dely*expfac;
	atoms->Fs(2,aj) -= delz*expfac;
	
	// need reverse comm here to get other contributions to integrand. should be fine
	
	}

    }

  }
  return;
}


/* set global cutoff  and default parameters */
void PairGridAtomGaussian::settings(std::vector<std::string> params)
{

  double cut = std::stod(params[0]);

  cutsq.resize(atoms->ntypes);
  epsilonstrength.resize(atoms->ntypes);
  nucwidth.resize(atoms->ntypes);
  phi.resize(atoms->ntypes);


  // set some default values
  for (int i = 0; i < epsilonstrength.size(); i++) {
    cutsq[i] = cut*cut;
    epsilonstrength[i] = 1.0;
    nucwidth[i] = 1.0;
    phi[i] = 0.1;
  }
  
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairGridAtomGaussian::coeff(std::vector<std::string> params)
{

  int type = std::stoi(params[0]);
  
  if (type > epsilonstrength.size() || type < 0)
    throw std::runtime_error("Pair grid/atom/gaussian atom type out of range.");
  
  

  params.erase(params.begin());

  if (params.size() > 6) 
    throw std::runtime_error("Pair grid/atom/gaussian invalid coeffs.");  
  
  int iarg = 0;
  while (iarg < params.size()) {

    if (params[iarg] == "epsilon") {
      epsilonstrength[type] = std::stod(params[iarg+1]);
      iarg += 2;
    } else if (params[iarg] == "nucwidth") {
      nucwidth[type] = std::stod(params[iarg+1]);
      iarg += 2;
    } else if (params[iarg] == "phi") {
      phi[type] = std::stod(params[iarg+1]);
      iarg += 2;
    } else
      throw std::runtime_error("Pair grid/atom/gaussian invalid coeffs.");
    
  }
  

}
