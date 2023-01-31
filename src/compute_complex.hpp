#ifndef PHAFD_COMPUTE_COMPLEX_HPP
#define PHAFD_COMPUTE_COMPLEX_HPP

#include <complex>
#include "compute.hpp"


namespace psPDE {
  template<typename>
  class fftw_MPI_3Darray;
}

namespace PHAFD_NS {

class ComputeComplex : public Compute
{
public:
  ComputeComplex(PHAFD *);


  virtual void init(const std::vector<std::string> &) override;
  virtual void in_fourier() override;
  virtual void end_of_step() override {};

  
private:
  psPDE::fftw_MPI_3Darray<std::complex<double>> *fftw3_arr;

  std::string which_quant;

  double prefac;
};
  

}
#endif
