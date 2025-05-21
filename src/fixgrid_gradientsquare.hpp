
#ifndef PHAFD_FIXGRID_GRADIENTSQUARE_HPP
#define PHAFD_FIXGRID_GRADIENTSQUARE_HPP


#include "fixgrid_laplacianphi.hpp"

#include <memory>



namespace PHAFD_NS {


class FixGridGradientSquare : public FixGridLaplacianPhi {
public:
  FixGridGradientSquare(PHAFD *);

  virtual void init(const std::vector<std::string> &) override;
  
  virtual void post_force() override;
  
private:


  double prefac;
};

}

#endif
