#ifndef PHAFD_COMPUTE_COMPLEX_HPP
#define PHAFD_COMPUTE_COMPLEX_HPP

#include <complex>
#include "compute.hpp"


namespace fftwArr {
  template<typename>
  class array3D;
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
  fftwArr::array3D<std::complex<double>> *fftw3_arr;

  std::string which_quant;

  double prefac;
};
  

}
#endif
