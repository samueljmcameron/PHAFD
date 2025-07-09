#include <algorithm>
#include <iostream>
#include "utility.hpp"

#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "integrate.hpp"
#include "fixgrid_divergence.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;


FixGridDivergence::FixGridDivergence(PHAFD *phafd) : Fix(phafd),
						     plan_set(false) {
  
};



void FixGridDivergence::init(const std::vector<std::string> &v_line,
			     bool add_to_names)
{

  Fix::init(v_line,add_to_names);

  if (v_line.size() > 1)
    throw std::runtime_error("invalid fix/divergence");

  divergence = nullptr;
  ft_divergence = nullptr;

  
}

FixGridDivergence::~FixGridDivergence()
{

  if (plan_set) {
    fftw_destroy_plan(backward_divergence);
  }
  
}

void FixGridDivergence::setup()
{


  normalization = (grid->boxgrid[0]*grid->boxgrid[1]*grid->boxgrid[2]);

  int Nx = grid->boxgrid[0];
  int Ny = grid->boxgrid[1];
  int Nz = grid->boxgrid[2];

  if (!divergence)
    divergence = std::make_unique<fftwArr::array3D<double>>
      (world,"divergence" +  name,Nx,Ny,Nz);
  if (!ft_divergence)
    ft_divergence
      = std::make_unique<fftwArr::array3D<std::complex<double>>>
      (world,"ft_divergence" +  name,Nx,Ny,Nz);
  
  backward_divergence
    = fftw_mpi_plan_dft_c2r_3d(Nz,Ny,Nx,
			       reinterpret_cast<fftw_complex*>
			       (ft_divergence->data()),
			       divergence->data(),world,
			       FFTW_MPI_TRANSPOSED_IN);
  
  

  plan_set = true;  
  
}


void FixGridDivergence::calculate_divergence(const
					     fftwArr::array3D<
					     std::complex<double>>
					     *fftw3_arr_x,
					     const
					     fftwArr::array3D<
					     std::complex<double>>
					     *fftw3_arr_y,
					     const
					     fftwArr::array3D<
					     std::complex<double>>
					     *fftw3_arr_z,
					     bool invert_fftw)

// std::array<*fftwArr::array3D<std::complex<double>>,SIZE> fftw3_arr
// or
// std::vector<*fftwArr::array3D<std::complex<double>>> fftw3_arr
{
  
  local0start = ft_divergence->get_local0start();
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
  
  for (int i = 0; i < ft_divergence->Nz(); i++) {
    
    if (i + local0start > globalNz/2) 
      l = -globalNz + i + local0start;
    else
      l = i + local0start;
    
    for (int j = 0; j < ft_divergence->Ny(); j++) {
      
      if (j > globalNy/2)
	m = -globalNy + j;
      else
	m = j;
      
      for (int k = 0; k < ft_divergence->Nx(); k++) {
	
	n = k;
	// due to transposes,
	// real space z is represented by qy
	// real space y is represented by qz

	(*ft_divergence)(i,j,k)
	  = (*fftw3_arr_x)(i,j,k)*idqx*n
	  + (*fftw3_arr_y)(i,j,k)*idqz*l
	  + (*fftw3_arr_z)(i,j,k)*idqy*m;
	
      }
    }
  }

  if (invert_fftw)
    fftw_execute(backward_divergence);
  
}


void FixGridDivergence::calculate_divergence(const
  std::array<std::unique_ptr<fftwArr::array3D<std::complex<double>>>,
					     3> &fftw3_arr,
					     bool invert_fftw)

// std::array<*fftwArr::array3D<std::complex<double>>,SIZE> fftw3_arr
// or
// std::vector<*fftwArr::array3D<std::complex<double>>> fftw3_arr
{
  
  local0start = ft_divergence->get_local0start();
  const int globalNy = grid->ft_boxgrid[1];
  const int globalNz = grid->ft_boxgrid[2];

  if (fftw3_arr.size() != 3)
    throw std::runtime_error("Error: calculating divergence on d!=3 object.");

  std::complex<double> idqx(0,domain->dqx());
  std::complex<double> idqy(0,domain->dqy());
  std::complex<double> idqz(0,domain->dqz());
  

  if (invert_fftw) {
    idqx /= normalization;
    idqy /= normalization;
    idqz /= normalization;
  }
  
  double l,m,n;
  
  for (int i = 0; i < ft_divergence->Nz(); i++) {
    
    if (i + local0start > globalNz/2) 
      l = -globalNz + i + local0start;
    else
      l = i + local0start;
    
    for (int j = 0; j < ft_divergence->Ny(); j++) {
      
      if (j > globalNy/2)
	m = -globalNy + j;
      else
	m = j;
      
      for (int k = 0; k < ft_divergence->Nx(); k++) {
	
	n = k;
	// due to transposes,
	// real space z is represented by qy
	// real space y is represented by qz

	(*ft_divergence)(i,j,k)
	  = (*fftw3_arr[0])(i,j,k)*idqx*n
	  + (*fftw3_arr[1])(i,j,k)*idqz*l
	  + (*fftw3_arr[2])(i,j,k)*idqy*m;
	
      }
    }
  }

  if (invert_fftw)
    fftw_execute(backward_divergence);
  
}

