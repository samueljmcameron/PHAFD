
#ifndef PHAFD_FIXGRID_VELOCITY_HPP
#define PHAFD_FIXGRID_VELOCITY_HPP


#include "fix.hpp"

#include <memory>

namespace fftwArr {
  template <typename T>
  class array3D;
}

namespace PHAFD_NS {

class FixGridVdet;
class FixGridVtherm;


class FixGridVelocity : public Fix {
public:
  FixGridVelocity(PHAFD *);

  virtual void init(const std::vector<std::string> &,bool) override;
  
  virtual void setup() override;
  virtual void start_of_step() override;  

  virtual void pre_final_integrate() override;
  
  virtual void post_final_integrate(bool /* flag */) override {};
  virtual void initial_integrate() override {};
  virtual void post_force() override {};

  virtual void final_integrate() override {};

  virtual void reset_dt() override;

  virtual void end_of_step() override {};
  
private:

  std::unique_ptr<FixGridVdet> vdet;
  std::unique_ptr<FixGridVtherm> vtherm;  
  

};

}

#endif
