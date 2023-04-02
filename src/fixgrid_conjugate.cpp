#include <algorithm>

#include "utility.hpp"


#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "integrate.hpp"
#include "fixgrid_conjugate.hpp"
#include "conjugate_volfrac.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;

template <typename T>
FixGridConjugate<T>::FixGridConjugate(PHAFD *phafd) : Fix(phafd) {};


template <typename T>
void FixGridConjugate<T>::init(const std::vector<std::string> &v_line)
{

  Fix::init(v_line);

  std::vector<std::string> new_v_line(v_line);


  int seed = std::stoi(new_v_line[1]);

  seed = utility::make_unique_seed(seed,world,commbrick->me,commbrick->nprocs);
  
  new_v_line.erase(new_v_line.begin(),new_v_line.begin()+2);

  new_v_line.push_back("seed");
  new_v_line.push_back(std::to_string(seed));


  conjugate = std::make_unique<T>(phafd);

  
  conjugate->readCoeffs(new_v_line);
  

}

template <typename T>  
void FixGridConjugate<T>::setup()
{

}
template <typename T>
void FixGridConjugate<T>::reset_dt()
{

  dt = integrate->dt;

  conjugate->reset_dt(dt);
  
}

template <typename T>
void FixGridConjugate<T>::start_of_step()
{

  // calculate the gradient of phi.
  fftw_execute(grid->forward_phi);

  const int local0start = grid->ft_phi->get_local0start();
  const int globalNy = grid->ft_boxgrid[1];
  const int globalNz = grid->ft_boxgrid[2];


  std::complex<double> idqx(0,domain->dqx());
  std::complex<double> idqy(0,domain->dqy());
  std::complex<double> idqz(0,domain->dqz());

  double l,m,n;

  for (int i = 0; i < grid->ft_phi->Nz(); i++) {
    
    if (i + local0start > globalNz/2) 
      l = -globalNz + i + local0start;
    else
      l = i + local0start;
    
    for (int j = 0; j < grid->ft_phi->Ny(); j++) {

      if (j > globalNy/2)
	m = -globalNy + j;
      else
	m = j;
      
      for (int k = 0; k < grid->ft_phi->Nx(); k++) {

	n = k;
	
	(*grid->ft_gradphi[0])(i,j,k) = (*grid->ft_phi)(i,j,k)*idqx*n;
	(*grid->ft_gradphi[1])(i,j,k) = (*grid->ft_phi)(i,j,k)*idqy*m;
	(*grid->ft_gradphi[2])(i,j,k) = (*grid->ft_phi)(i,j,k)*idqz*l;

	
      }
    }
  }

  fftw_execute(grid->backward_gradphi[0]);
  fftw_execute(grid->backward_gradphi[1]);
  fftw_execute(grid->backward_gradphi[2]);

  double factor = grid->boxgrid[0]*grid->boxgrid[1]*grid->boxgrid[2];
  (*grid->gradphi[0]) /= factor;
  (*grid->gradphi[1]) /= factor;
  (*grid->gradphi[2]) /= factor;
  
}


template <typename T>
void FixGridConjugate<T>::pre_final_integrate()
{


  fftw_execute(grid->forward_chempot);

  didnotintegrate = true;
  return;
}

template <typename T>
void FixGridConjugate<T>::final_integrate()
{


  conjugate->update();
  didnotintegrate = false;
  
  return;
}

template <typename T>
void FixGridConjugate<T>::post_final_integrate()
{

  fftw_execute(grid->backward_phi);

  if (didnotintegrate) {
    double factor = grid->boxgrid[0]*grid->boxgrid[1]*grid->boxgrid[2];
    
    (*grid->phi) /= factor;
  }
  
  return;
}


template class FixGridConjugate<ConjugateVolFrac>;
