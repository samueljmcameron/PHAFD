#include <algorithm>

#include "utility.hpp"

#include "ps_pde/conjugate_volfrac.hpp"

#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"

#include "fixgrid_conjugate.hpp"

using namespace PHAFD_NS;

template <typename T>
FixGridConjugate<T>::FixGridConjugate(PHAFD *phafd) : FixGrid(phafd) {};


template <typename T>
void FixGridConjugate<T>::init(const std::vector<std::string> &v_line)
{

  FixGrid::init(v_line);

  std::vector<std::string> new_v_line(v_line);


  int seed = std::stoi(new_v_line[1]);

  seed = utility::make_unique_seed(seed,world,commbrick->me,commbrick->nprocs);
  
  new_v_line.erase(new_v_line.begin(),new_v_line.begin()+2);

  new_v_line.push_back("seed");
  new_v_line.push_back(std::to_string(seed));


  conjugate = std::make_unique<T>(*(domain->ps_domain),*(grid->ps_grid));


  
  conjugate->readCoeffs(new_v_line);
  

}

template <typename T>  
void FixGridConjugate<T>::setup()
{

}
template <typename T>
void FixGridConjugate<T>::reset_dt(double timestep)
{

  dt = timestep;

  conjugate->reset_dt(dt);
  
}

template <typename T>  
void FixGridConjugate<T>::initial_integrate()
{
  
}
template <typename T>
void FixGridConjugate<T>::final_integrate()
{

  fftw_execute(grid->ps_grid->forward_phi);
  fftw_execute(grid->ps_grid->forward_nonlinear);

  conjugate->update();

  fftw_execute(grid->ps_grid->backward_phi);
  
  return;
}


template class FixGridConjugate<psPDE::ConjugateVolFrac>;
