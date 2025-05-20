
#ifndef PHAFD_FIXGRID_V_DOT_GRADPHI_HPP
#define PHAFD_FIXGRID_V_DOT_GRADPHI_HPP


#include "fix.hpp"

#include <memory>

namespace fftwArr {
  template <typename T>
  class array3D;
}

namespace PHAFD_NS {

class FixGridVelocity;



class FixGridVdotGradPhi : public Fix {
public:
  FixGridVdotGradPhi(PHAFD *);

  virtual void init(const std::vector<std::string> &) override;
  
  virtual void setup() override;
  virtual void start_of_step() override;  

  virtual void pre_final_integrate() override;
  virtual void post_final_integrate() override;
  virtual void initial_integrate() override {};
  virtual void post_force() override {};

  virtual void final_integrate() override {};

  virtual void reset_dt() override;

  virtual void end_of_step() override {};
  
private:

  std::unique_ptr<FixGridVelocity> velocity;
  

};

}

#endif
