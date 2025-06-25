
#ifndef PHAFD_FIXGRID_DIVERGENCE_HPP
#define PHAFD_FIXGRID_DIVERGENCE_HPP


#include "fix.hpp"

#include <memory>


namespace fftwArr {
  template<typename>
  class array3D;
}

namespace PHAFD_NS {

class FixGridDivergence : public Fix {
public:
  FixGridDivergence(PHAFD *);
  ~FixGridDivergence();

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

  //template <class arr_type>
  void calculate_divergence(const std::array<std::unique_ptr<
			    fftwArr::array3D<std::complex<double>>>,3> &,
			    bool invert_fftw=true);
  std::unique_ptr<fftwArr::array3D<std::complex<double>>> ft_divergence;
  std::unique_ptr<fftwArr::array3D<double>> divergence;
  
private:

  double normalization;
  bool plan_set;

  


  fftw_plan backward_divergence;

  
};

}

#endif
