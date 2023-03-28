#include "compute_complex.hpp"

#include <cmath>

#include "grid.hpp"
#include "ps_pde/fftw_mpi_3darray.hpp"

using namespace PHAFD_NS;


ComputeComplex::ComputeComplex(PHAFD *phafd) : Compute(phafd) {
  per_ftgrid = true;
}


void ComputeComplex::init(const std::vector<std::string> &v_line) {

  Compute::init(v_line);


  std::string arrname = v_line.at(1);

  if (arrname == "ft_phi") {
    fftw3_arr = grid->ft_phi;
  } else if (arrname == "ft_chempot") {
    fftw3_arr = grid->ft_chempot;
  } else if (arrname == "ft_noise") {
    fftw3_arr = grid->ft_noise;
  } else {
    throw std::runtime_error("Invalid array to compute in compute modulus.");
  }

  if (fftw3_arr == nullptr)
    throw std::runtime_error("Cannot compute complex of array " + arrname
			     + std::string(", not allocated."));


  
  which_quant = v_line.at(2);

  if (which_quant != "real" && which_quant != "imag" && which_quant != "modulus"
      && which_quant != "norm")
    throw std::runtime_error("Invalid argument in compute complex.");

  
  Nz = fftw3_arr->Nz();
  Ny = fftw3_arr->Ny();
  Nx = fftw3_arr->Nx();

  prefac = 1./(grid->dx()*grid->dy()*grid->dz());
  
  array.resize(Nz*Ny*Nx);
  
}

void ComputeComplex::in_fourier()
{

  if (!this_step) return;
  int count = 0;

  if (which_quant == "modulus") {
  
    for (int i = 0; i < Nz; i++)
      for (int j = 0; j < Ny; j++)
	for (int k = 0; k < Nx; k++)
	  array[count++] = std::abs((*fftw3_arr)(i,j,k))*prefac;

  } else if (which_quant == "norm") {
  
    for (int i = 0; i < Nz; i++)
      for (int j = 0; j < Ny; j++)
	for (int k = 0; k < Nx; k++) 
	  array[count++] = std::norm((*fftw3_arr)(i,j,k))*prefac*prefac;


  } else if (which_quant == "real") {
  
    for (int i = 0; i < Nz; i++)
      for (int j = 0; j < Ny; j++)
	for (int k = 0; k < Nx; k++)
	  array[count++] = (*fftw3_arr)(i,j,k).real()*prefac;

  } else if (which_quant == "imag") {
  
    for (int i = 0; i < Nz; i++)
      for (int j = 0; j < Ny; j++)
	for (int k = 0; k < Nx; k++)
	  array[count++] = (*fftw3_arr)(i,j,k).imag()*prefac;

  }
  
}
