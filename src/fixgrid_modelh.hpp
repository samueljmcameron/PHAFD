
#ifndef PHAFD_FIXGRID_MODELH_HPP
#define PHAFD_FIXGRID_MODELH_HPP


#include "fix.hpp"

#include <memory>

namespace fftwArr {
  template<typename>
  class array3D;
}

namespace PHAFD_NS {

class ConjugateNoise;


class FixGridModelH : public Fix {
public:
  FixGridModelH(PHAFD *);

  virtual void init(const std::vector<std::string> &) override;
  
  virtual void setup() override;
  virtual void reset_dt() override;
  
  virtual void start_of_step() override;
  virtual void initial_integrate() override;
  virtual void post_force() override ;
  virtual void pre_final_integrate() override;
  virtual void final_integrate() override;
  virtual void post_final_integrate() override;
  virtual void end_of_step() override ;

  
private:

  std::unique_ptr<ConjugateNoise> conjugate;

  std::vector<std::unique_ptr<Fix>> local_fixes; // store velocity and gradphi fixes.

  double normalization;
  double mobility;

  void point_update(int,int,int);

  void origin_update();


  
  
};

}

#endif
