
#ifndef PHAFD_FIXGRID_LAPLACIAN_HPP
#define PHAFD_FIXGRID_LAPLACIAN_HPP


#include "fix.hpp"

#include <memory>


namespace fftwArr {
  template<typename>
  class array3D;
}

namespace PHAFD_NS {

class FixGridLaplacian : public Fix {
public:
  FixGridLaplacian(PHAFD *);
  ~FixGridLaplacian();

  virtual void init(const std::vector<std::string> &) override;

  virtual void setup() override;
  virtual void pre_final_integrate() override {};
  virtual void post_final_integrate() override {};
  virtual void start_of_step() override {};
  virtual void initial_integrate() override {};
  virtual void post_force() override {};
  virtual void final_integrate() override {};
  virtual void end_of_step() override {};

  //

  void calculate_laplacian(const fftwArr::array3D<std::complex<double>> *);
  std::unique_ptr<fftwArr::array3D<double>> laplacian;
  
private:

  bool plan_set;

  std::unique_ptr<fftwArr::array3D<std::complex<double>>> ft_laplacian;
  
  fftw_plan backward_laplacian;

  
  
  
};

}

#endif
