
#ifndef PHAFD_FIXGRID_MODELB_ADVECT_HPP
#define PHAFD_FIXGRID_MODELB_ADVECT_HPP


#include "fixgrid_modelb.hpp"

#include <memory>

namespace fftwArr {
  template<typename>
  class array3D;
}

namespace PHAFD_NS {
class FixGridGradient;
class ConjugateNoise;


class FixGridModelBadvect : public FixGridModelB {
public:
  FixGridModelBadvect(PHAFD *);

  ~FixGridModelBadvect();
  virtual void init(const std::vector<std::string> &,bool) override;
  
  virtual void setup() override;

  virtual void start_of_step() override;
  
  virtual void post_final_integrate(bool invert_fft=true) override;
  
private:
  
  std::unique_ptr<FixGridGradient> gradphifix;

  std::array<std::unique_ptr<fftwArr::array3D<double>>,3> flow;

  void set_flow();


  
  
};

}

#endif
