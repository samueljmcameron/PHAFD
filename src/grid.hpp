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
  psPDE::fftw_MPI_3Darray<double> *chempot;
  psPDE::fftw_MPI_3Darray<double> *gradphi_x,*gradphi_y,*gradphi_z;
  

  psPDE::fftw_MPI_3Darray<std::complex<double>> *ft_phi;
  psPDE::fftw_MPI_3Darray<std::complex<double>> *ft_chempot;
  psPDE::fftw_MPI_3Darray<std::complex<double>> *ft_noise;
  psPDE::fftw_MPI_3Darray<std::complex<double>> *ft_gradphi_x,*ft_gradphi_y,*ft_gradphi_z;
  
  std::unique_ptr<psPDE::Grid> ps_grid;  
  bool gridset,gridpopulated;


  Grid(PHAFD *);
  ~Grid();
  void create(const std::vector<std::string> &);
  void populate(const std::vector<std::string> &);

  double dx() const;
  double dy() const;
  double dz() const;
  
private:


  void sinusoidal_grid(const std::string &,double);


  

};

}


#endif
