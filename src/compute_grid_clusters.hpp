#ifndef PHAFD_COMPUTE_GRID_CLUSTERS_HPP
#define PHAFD_COMPUTE_GRID_CLUSTERS_HPP


#include "compute.hpp"

namespace fftwArr {

  template<typename>
  class array3D;

}

namespace PHAFD_NS {

class ComputeGridClusters : public Compute
{
public:
  ComputeGridClusters(PHAFD *);


  virtual void init(const std::vector<std::string> &) override;

  
  virtual void in_fourier() override {};
  virtual void end_of_step() override;

  
private:
  fftwArr::array3D<double> *fftw3_arr;

  std::string condition;
  std::vector<int> intarray,lefts,rights;

  double threshold;


  void send_to_neighbors();
  
};
  

}
#endif
