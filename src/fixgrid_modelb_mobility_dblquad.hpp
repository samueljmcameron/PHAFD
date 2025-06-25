
#ifndef PHAFD_FIXGRID_MODELB_MOBILITY_DBLQUAD_HPP
#define PHAFD_FIXGRID_MODELB_MOBILITY_DBLQUAD_HPP


#include "fixgrid_modelb_mobility_base.hpp"


namespace fftwArr {
  template<typename>
  class array3D;
}

namespace PHAFD_NS {

class FixGridModelBMobilityDBLquad : public FixGridModelBMobilityBase {
public:
  FixGridModelBMobilityDBLquad(PHAFD *);

  virtual void init(const std::vector<std::string> &) override;
  
  virtual void calculate_sqrt_mobility() override;

private:
  double prefactor;

  
};

}

#endif
