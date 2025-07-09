
#ifndef PHAFD_FIXGRID_SPHERICAL_AVERAGE_HPP
#define PHAFD_FIXGRID_SPHERICAL_AVERAGE_HPP


#include "fix.hpp"
#include <fftw3-mpi.h>

#include <memory>
#include <complex>

namespace fftwArr {
  template<typename>
  class array3D;
}

namespace PHAFD_NS {


class FixGridSphericalAverage : public Fix {
public:
  FixGridSphericalAverage(PHAFD *);


  virtual void init(const std::vector<std::string> &) override;
  
  virtual void setup() override;

  virtual void start_of_step() override;
  
  virtual void initial_integrate() override {};
  virtual void post_force() override {};
  virtual void pre_final_integrate() override {};
  virtual void final_integrate() override {};
  virtual void post_final_integrate() override {};
  
  virtual void reset_dt() override {};

  virtual void end_of_step() override;
  
private:

  std::vector<int> counts;
  std::vector<int> global_counts;
  std::vector<double> output;

  void FixGridSphericalAverage::loop(const fftwArr::array3D<
				     std::complex<double>> *);

  
};

}

#endif
