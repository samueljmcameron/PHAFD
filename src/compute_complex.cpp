#include "compute_complex.hpp"

#include <cmath>
#include <iostream>
#include "utility.hpp"
#include "grid.hpp"
#include "fix.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;


ComputeComplex::ComputeComplex(PHAFD *phafd) : Compute(phafd) {
  per_ftgrid = true;
}

ComputeComplex::~ComputeComplex()
{
}

void ComputeComplex::init(const std::vector<std::string> &v_line) {

  Compute::init(v_line);


  std::string arrname = v_line.at(1);

  if (arrname == "ft_phi") {
    fftw3_arr = grid->ft_phi.get();
  } else if (arrname == "ft_chempot") {
    fftw3_arr = grid->ft_chempot.get();
  } else if (arrname == "ft_vtherm_x") {
    fftw3_arr = grid->ft_vtherm[0].get();
  } else if (arrname == "ft_vtherm_y") {
    fftw3_arr = grid->ft_vtherm[1].get();
  } else if (arrname == "ft_vtherm_z") {
    fftw3_arr = grid->ft_vtherm[2].get();
  } else if (arrname == "ft_vdet_x") {
    fftw3_arr = grid->ft_vdet[0].get();
  } else if (arrname == "ft_vdet_y") {
    fftw3_arr = grid->ft_vdet[1].get();
  } else if (arrname == "ft_vdet_z") {
    fftw3_arr = grid->ft_vdet[2].get();
  } else if (arrname == "ft_v_dot_gradphi") {
    fftw3_arr = grid->ft_v_dot_gradphi.get();
  }  else {
    utility::find_array_component(arrname,phafd,fftw3_arr);
  }
  if (fftw3_arr == nullptr)
    throw std::runtime_error("Cannot compute complex of array " + arrname
			     + std::string(", not allocated."));

  
  which_quant = v_line.at(2);

  if (which_quant != "real" && which_quant != "imag" && which_quant != "modulus"
      && which_quant != "norm")
    throw std::runtime_error("Invalid argument in compute complex.");


  localNz = fftw3_arr->Nz();
  localNy = fftw3_arr->Ny();
  localNx = fftw3_arr->Nx();
  prefac = 1.0;
  // too much confusion using prefactor that isn't unity.


  if (!output) 
    output = std::make_unique<fftwArr::array3D<double>
			      >(world,name + std::string("_output"),
				localNx,
				grid->boxgrid[1],
				grid->boxgrid[2]);

  array.resize(localNx*localNy*localNz);
  realFFTWarray.push_back(output.get());

}


void ComputeComplex::in_fourier()
{

  if (!this_step) return;
  int count = 0;
  if (which_quant == "modulus") {
  
    for (int i = 0; i < localNz; i++)
      for (int j = 0; j < localNy; j++)
	for (int k = 0; k < localNx; k++)
	  (*output)(i,j,k) = std::abs((*fftw3_arr)(i,j,k))*prefac;
	  //array[count++] = std::abs((*fftw3_arr)(i,j,k))*prefac;

  } else if (which_quant == "norm") {
  
    for (int i = 0; i < localNz; i++)
      for (int j = 0; j < localNy; j++)
	for (int k = 0; k < localNx; k++) 
	  (*output)(i,j,k) = std::norm((*fftw3_arr)(i,j,k))*prefac*prefac;
    //array[count++] = std::norm((*fftw3_arr)(i,j,k))*prefac*prefac;

  } else if (which_quant == "real") {
  
    for (int i = 0; i < localNz; i++)
      for (int j = 0; j < localNy; j++)
	for (int k = 0; k < localNx; k++)
	  (*output)(i,j,k) = (*fftw3_arr)(i,j,k).real()*prefac;
    //array[count++] = (*fftw3_arr)(i,j,k).real()*prefac;
  } else if (which_quant == "imag") {
  
    for (int i = 0; i < localNz; i++)
      for (int j = 0; j < localNy; j++)
	for (int k = 0; k < localNx; k++)
	  (*output)(i,j,k) = (*fftw3_arr)(i,j,k).imag()*prefac;
    //array[count++] = (*fftw3_arr)(i,j,k).imag()*prefac;
  }
  
}
