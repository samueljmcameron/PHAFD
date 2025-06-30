
#ifndef PHAFD_FIXGRID_MODELB_MOBILITY_CONSTANT_HPP
#define PHAFD_FIXGRID_MODELB_MOBILITY_CONSTANT_HPP


#include "fixgrid_modelb_mobility_base.hpp"


namespace fftwArr {
  template<typename>
  class array3D;
}

namespace PHAFD_NS {

class FixGridModelBMobilityConstant : public FixGridModelBMobilityBase {
public:
  FixGridModelBMobilityConstant(PHAFD *);

  virtual void init(const std::vector<std::string> &) override;
  
  virtual void calculate_sqrt_mobility() override;
  virtual void calculate_mobility_deriv() override;

private:
  double prefactor;

  
};

}

#endif
