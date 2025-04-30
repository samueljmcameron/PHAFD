
#ifndef PHAFD_COMPUTE_GRID_CHEMPOTGRADPHI_HPP
#define PHAFD_COMPUTE_GRID_CHEMPOTGRADPHI_HPP


#include "compute.hpp"

#include <memory>

namespace fftwArr {
  template <typename T>
  class array3D;
}

namespace PHAFD_NS {

class ConjugateNoise;

class ComputeGridChemPotGradPhi : public Compute {
public:
  ComputeGridChemPotGradPhi(PHAFD *);

  virtual void init(const std::vector<std::string> &) override;
  
  virtual void setup() override;

  virtual void start_of_step() override;
  
  virtual void initial_integrate() override {};
  virtual void post_force() override {};
  virtual void pre_final_integrate() override {};
  virtual void final_integrate() override {};
  virtual void post_final_integrate() override {}  
  virtual void reset_dt() override;

  virtual void end_of_step() override {};
  
private:

  void compute_vtherm();
  void set_vtherm(int, int, int) ;

  std::array<fftwArr::array3D<std::complex<double>>,3> chempotgradphi;

  
};

}

#endif
