
#ifndef PHAFD_FIXGRID_GRADIENTSQUARE_HPP
#define PHAFD_FIXGRID_GRADIENTSQUARE_HPP


#include "fix.hpp"

#include <memory>



namespace PHAFD_NS {

class FixGridLaplacian;
  
class FixGridGradientSquare : public Fix {
public:
  FixGridGradientSquare(PHAFD *);

  virtual void init(const std::vector<std::string> &) override;
  
  virtual void post_force() override;


  virtual void setup() override;
  virtual void start_of_step() override {};

  virtual void pre_final_integrate() override {};
  virtual void post_final_integrate(bool /* flag */) override {};
  virtual void initial_integrate() override {};


  virtual void final_integrate() override {};

  virtual void reset_dt() override {};

  virtual void end_of_step() override {};

  
  
private:

  std::unique_ptr<FixGridLaplacian> laplacian;

  double prefac;
  bool ft_phi_flag;
};
  
}

#endif
