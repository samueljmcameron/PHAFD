
#ifndef PHAFD_FIXGRID_GRADPHI_HPP
#define PHAFD_FIXGRID_GRADPHI_HPP


#include "fix.hpp"

#include <memory>

namespace PHAFD_NS {

class FixGridGradPhi : public Fix {
public:
  FixGridGradPhi(PHAFD *);

  virtual void init(const std::vector<std::string> &) override;
  
  virtual void setup() override;

  virtual void start_of_step() override;
  
  virtual void initial_integrate() override {};
  virtual void post_force() override {};
  virtual void pre_final_integrate() override {};
  virtual void final_integrate() override {};
  virtual void post_final_integrate() override {};
  

  virtual void end_of_step() override {};
private:


  bool once;
};

}

#endif
