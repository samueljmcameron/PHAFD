#ifndef PHAFD_COMPUTE_FT_SPHERICAL_BIN_HPP
#define PHAFD_COMPUTE_FT_SPHERICAL_BIN_HPP

#include <vector>
#include "compute.hpp"


namespace fftwArr {
  template<typename>
  class array3D;
}

namespace PHAFD_NS {

class Fix;
  
class ComputeFtSphericalBin : public Compute
{
public:
  ComputeFtSphericalBin(PHAFD *);
  ~ComputeFtSphericalBin();


  virtual void init(const std::vector<std::string> &,bool) override;
  virtual void in_fourier() override {};
  virtual void start_of_step() override;
  virtual void end_of_step() override;

  
  
private:

  template <int Tp_COUNT> void loop();

  Compute *compute;
  Fix *fix;

  double *input_array;
  
  int input_component;   // which component of input_array is being used
  int input_nc; // number of components for input array

  std::vector<double> output;
  std::vector<int> counts, global_counts;

  int nbins;
  double dqbins;
};
  

}
#endif
