
#ifndef PHAFD_FIXGRID_GRADIENT_HPP
#define PHAFD_FIXGRID_GRADIENT_HPP


#include "fix.hpp"

#include <memory>


namespace fftwArr {
  template<typename>
  class array3D;
}

namespace PHAFD_NS {

class FixGridGradient : public Fix {
public:
  FixGridGradient(PHAFD *);
  ~FixGridGradient();

  virtual void init(const std::vector<std::string> &) override;

  virtual void setup() override;
  virtual void pre_final_integrate() override {};
  virtual void post_final_integrate() override {};
  virtual void start_of_step() override {};
  virtual void initial_integrate() override;
  virtual void post_force() override {};
  virtual void final_integrate() override {};
  virtual void end_of_step() override {};

  //

  void calculate_gradient(const fftwArr::array3D<std::complex<double>> *,
			  bool invert_fftw=true);
  std::array<std::unique_ptr<fftwArr::array3D<double>>,3> gradient;
  std::array<std::unique_ptr<fftwArr::array3D<std::complex<double>>>,3>
  ft_gradient;
  
private:

  bool plan_set;

  double normalization;

  


  std::array<fftw_plan,3> backward_gradient;

  
  
  
};

}

#endif
