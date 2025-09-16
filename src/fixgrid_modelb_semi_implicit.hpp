
#ifndef PHAFD_FIXGRID_MODELB_SEMI_IMPLICIT_HPP
#define PHAFD_FIXGRID_MODELB_SEMI_IMPLICIT_HPP


#include "fix.hpp"

#include <memory>

namespace fftwArr {
  template<typename>
  class array3D;
}

namespace PHAFD_NS {
class FixGridGradient;
class ConjugateNoise;


class FixGridModelBSemiImplicit : public Fix {
public:
  FixGridModelBSemiImplicit(PHAFD *);

  virtual void init(const std::vector<std::string> &,bool) override;
  
  virtual void setup() override;

  virtual void start_of_step() override;
  
  virtual void initial_integrate() override {};
  virtual void post_force() override {};
  virtual void pre_final_integrate() override;
  virtual void final_integrate() override;
  virtual void post_final_integrate(bool invert_fft=true) override;
  
  virtual void reset_dt() override;

  virtual void end_of_step() override {};
private:
  
  std::unique_ptr<FixGridGradient> gradfix;

  std::unique_ptr<ConjugateNoise> conjugate;

  double normalization;
  double mobility;
  double kappa;

  void point_update(int,int,int);

  void origin_update();


  
  
};

}

#endif
