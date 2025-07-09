#include <algorithm>
#include <iostream>
#include "utility.hpp"

#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "integrate.hpp"
#include "fixgrid_laplacian.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;


FixGridLaplacian::FixGridLaplacian(PHAFD *phafd) : Fix(phafd),
						   plan_set(false) {
  
};



void FixGridLaplacian::init(const std::vector<std::string> &v_line,
			    bool add_to_names)
{

  Fix::init(v_line,add_to_names);

  if (v_line.size() > 1)
    throw std::runtime_error("invalid fix/laplacian");    

  laplacian = nullptr;
  ft_laplacian = nullptr;
  
}

FixGridLaplacian::~FixGridLaplacian()
{

  if (plan_set)
    fftw_destroy_plan(backward_laplacian);
    
}

void FixGridLaplacian::setup()
{

  
  int Nx = grid->boxgrid[0];
  int Ny = grid->boxgrid[1];
  int Nz = grid->boxgrid[2];

  if (!laplacian)
    laplacian = std::make_unique<fftwArr::array3D<double>>
      (world,"laplacian"+name,Nx,Nz,Ny);

  
  if (!ft_laplacian)
    ft_laplacian
      = std::make_unique<fftwArr::array3D<std::complex<double>>>
      (world,"ft_laplacian"+name,Nx,Nz,Ny);

  backward_laplacian = fftw_mpi_plan_dft_c2r_3d(Nz,Ny,Nx,
						reinterpret_cast<fftw_complex*>
						(ft_laplacian->data()),
						laplacian->data(),world,
						FFTW_MPI_TRANSPOSED_IN);

  plan_set = true;  
  
}

void FixGridLaplacian::calculate_laplacian(const fftwArr::array3D<std::complex<double>>
					   *fftw3_arr)
{


  double normalization
    = (grid->boxgrid[0]*grid->boxgrid[1]*grid->boxgrid[2]);

  
  local0start = fftw3_arr->get_local0start();
  const int globalNy = grid->ft_boxgrid[1];
  const int globalNz = grid->ft_boxgrid[2];
  
  double dqy = domain->dqy();
  double dqz = domain->dqz();
  double dqx = domain->dqx();
  
  dqx /= sqrt(normalization);
  dqy /= sqrt(normalization);
  dqz /= sqrt(normalization);
  
  double l,m,n;
  double neg_q2;
 
    
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
	
	neg_q2 = -(dqx*n*dqx*n + dqz*m*dqz*m + dqy*l*dqy*l);
	
	
	(*ft_laplacian)(i,j,k) = (*fftw3_arr)(i,j,k)*neg_q2;
	
      }
    }
  }
  
  fftw_execute(backward_laplacian);
  
}
