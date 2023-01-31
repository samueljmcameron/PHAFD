#ifndef PHAFD_GRID_HPP
#define PHAFD_GRID_HPP

#include <complex>
#include "pointers.hpp"

namespace psPDE {

  template<typename>
  class fftw_MPI_3Darray;
  class Grid;
}

namespace PHAFD_NS {


class Grid : protected Pointers
{
public:

  // points to psPDE::grid std::arrays (of size three)
  ptrdiff_t *boxgrid; // global number of grid points {Nx,Ny,Nz}
  ptrdiff_t *ft_boxgrid; // global number of grid points {Nx,Ny,Nz}

  psPDE::fftw_MPI_3Darray<double> *phi;
  psPDE::fftw_MPI_3Darray<double> *nonlinear;

  psPDE::fftw_MPI_3Darray<std::complex<double>> *ft_phi;
  psPDE::fftw_MPI_3Darray<std::complex<double>> *ft_nonlinear;
  psPDE::fftw_MPI_3Darray<std::complex<double>> *ft_noise;
  
  std::unique_ptr<psPDE::Grid> ps_grid;  



  Grid(PHAFD *);
  ~Grid();
  void create(const std::vector<std::string> &);
  void populate(const std::vector<std::string> &);

  double dx() const;
  double dy() const;
  double dz() const;
  
private:




  

};

}


#endif
