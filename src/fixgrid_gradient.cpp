#include <algorithm>
#include <iostream>
#include "utility.hpp"

#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "integrate.hpp"
#include "fixgrid_gradient.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;


FixGridGradient::FixGridGradient(PHAFD *phafd) : Fix(phafd),
						 plan_set(false) {
  
};



void FixGridGradient::init(const std::vector<std::string> &v_line,
			   bool add_to_names)
{

  Fix::init(v_line,add_to_names);

  if (v_line.size() > 1)
    throw std::runtime_error("invalid fix/gradient");

  for (int i = 0; i < 3; i++) {
    gradient[i] = nullptr;
    ft_gradient[i] = nullptr;
  }
  
}

FixGridGradient::~FixGridGradient()
{

  if (plan_set) {
    for (int i = 0; i < 3; i++)
      fftw_destroy_plan(backward_gradient[i]);
  }
  
}

void FixGridGradient::setup()
{

  normalization = (grid->boxgrid[0]*grid->boxgrid[1]*grid->boxgrid[2]);
  
  int Nx = grid->boxgrid[0];
  int Ny = grid->boxgrid[1];
  int Nz = grid->boxgrid[2];

  for (int i = 0; i < 3; i++) {
    if (!gradient[i])
      gradient[i] = std::make_unique<fftwArr::array3D<double>>
	(world,"gradient[" + std::to_string(i) + "]" +  name,Nx,Ny,Nz);
    if (!ft_gradient[i])
      ft_gradient[i]
	= std::make_unique<fftwArr::array3D<std::complex<double>>>
	(world,"ft_gradient[" + std::to_string(i) + "]" +  name,Nx,Ny,Nz);
    
    backward_gradient[i]
      = fftw_mpi_plan_dft_c2r_3d(Nz,Ny,Nx,
				 reinterpret_cast<fftw_complex*>
				 (ft_gradient[i]->data()),
				 gradient[i]->data(),world,
				 FFTW_MPI_TRANSPOSED_IN);
    
    
  }
  
  plan_set = true;  
  
}

void FixGridGradient::initial_integrate()
{
  calculate_gradient(grid->ft_phi.get());
  for (int dim = 0; dim < 3; dim ++)
    *grid->gradphi[dim] = *gradient[dim];
}

void FixGridGradient::calculate_gradient(const
					 fftwArr::array3D<
					 std::complex<double>>
					 *fftw3_arr,
					 bool invert_fftw) {
  
  local0start = fftw3_arr->get_local0start();
  const int globalNy = grid->ft_boxgrid[1];
  const int globalNz = grid->ft_boxgrid[2];
  
  
  std::complex<double> idqx(0,domain->dqx());
  std::complex<double> idqy(0,domain->dqy());
  std::complex<double> idqz(0,domain->dqz());
  

  if (invert_fftw) {
    idqx /= normalization;
    idqy /= normalization;
    idqz /= normalization;
  }
  
  double l,m,n;
  
  for (int i = 0; i < fftw3_arr->Nz(); i++) {
    
    if (i + local0start > globalNz/2) 
      l = -globalNz + i + local0start;
    else
      l = i + local0start;
    
    for (int j = 0; j < fftw3_arr->Ny(); j++) {
      
      if (j > globalNy/2)
	m = -globalNy + j;
      else
	m = j;
      
      for (int k = 0; k < fftw3_arr->Nx(); k++) {
	
	n = k;
	// due to transposes,
	// real space z is represented by qy
	// real space y is represented by qz
	(*ft_gradient[0])(i,j,k) = (*fftw3_arr)(i,j,k)*idqx*n;
	(*ft_gradient[1])(i,j,k) = (*fftw3_arr)(i,j,k)*idqz*l;
	(*ft_gradient[2])(i,j,k) = (*fftw3_arr)(i,j,k)*idqy*m;
	
      }
    }
  }

  if (invert_fftw)
    for (int i = 0; i < 3; i++) 
      fftw_execute(backward_gradient[i]);
  
}
