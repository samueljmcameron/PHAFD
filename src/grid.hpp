#ifndef PHAFD_GRID_HPP
#define PHAFD_GRID_HPP

#include <complex>
#include <array>
#include <fftw3-mpi.h>
#include "pointers.hpp"

namespace fftwArr {

  template<typename>
  class array3D;

}

namespace PHAFD_NS {


class Grid : protected Pointers
{
public:

  // points to psPDE::grid std::arrays (of size three)
  std::array<ptrdiff_t,3> boxgrid; // global number of grid points {Nx,Ny,Nz}
  std::array<ptrdiff_t,3> ft_boxgrid; // global number of grid points {Nx,Ny,Nz}


  /* concentration (model B) arrays */
  std::unique_ptr<fftwArr::array3D<double>> phi; // concentration field
  std::unique_ptr<fftwArr::array3D<double>> chempot; // chemical potential
  std::array<std::unique_ptr<fftwArr::array3D<double>>,3> gradphi; // gradients

  // and their fourier transforms
  std::unique_ptr<fftwArr::array3D<std::complex<double>>> ft_phi;
  std::unique_ptr<fftwArr::array3D<std::complex<double>>> ft_chempot;
  std::array<std::unique_ptr<fftwArr::array3D<std::complex<double>>>,3> ft_gradphi;

  // and their fftw plans
  fftw_plan forward_phi, backward_phi;
  fftw_plan forward_chempot, backward_chempot;
  std::array<fftw_plan,3> backward_gradphi;


  /* velocity (Navier-stokes) arrays */
  std::array<std::unique_ptr<fftwArr::array3D<double>>,3> velocity; // fluid velocity
  std::array<std::unique_ptr<fftwArr::array3D<double>>,3> vtherm; // fluid thermal noise
  
  // and their fourier transforms
  std::array<std::unique_ptr<fftwArr::array3D<std::complex<double>>>,3> ft_velocity;
  std::array<std::unique_ptr<fftwArr::array3D<std::complex<double>>>,3> ft_vtherm;

  std::array<fftw_plan,3> forward_velocity, backward_velocity;
  std::array<fftw_plan,3> forward_vtherm, backward_vtherm;
  

  /* Model H arrays (not already established) */

  std::unique_ptr<fftwArr::array3D<double>> vtherm_dot_gradphi;
  std::unique_ptr<fftwArr::array3D<std::complex<double>>> ft_vtherm_dot_gradphi;

  std::array<std::unique_ptr<fftwArr::array3D<std::complex<double>>>,3> ft_Znoise;

  
  bool gridset,gridpopulated;

  bool phi_set,velocity_set,modelH_set;


  Grid(PHAFD *);
  ~Grid();
  void create(const std::vector<std::string> &);
  void populate(const std::vector<std::string> &);

  double dx() const;
  double dy() const;
  double dz() const;
  
private:


  void sinusoidal_grid(const std::string &,double);
  void create_concentration(int , int , int);
  void create_noise(int, int, int);
  void noisy_constant(fftwArr::array3D<double>*, double,double, int);
  void constant_noise(double,double,int);
  void create_modelH(int, int, int);
  void flat_interface(fftwArr::array3D<double>*, double,double,double,double, int);
  void sphere(fftwArr::array3D<double>*,double, double , double ,
	      double);

  double sphericalshape(double , double, double);

  

};

}


#endif
