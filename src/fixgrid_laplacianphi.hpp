
#ifndef PHAFD_FIXGRID_LAPLACIANPHI_HPP
#define PHAFD_FIXGRID_LAPLACIANPHI_HPP


#include "fix.hpp"

#include <memory>

namespace PHAFD_NS {

class FixGridLaplacianPhi : public Fix {
public:
  FixGridLaplacianPhi(PHAFD *);

  virtual void init(const std::vector<std::string> &) override;
  
  virtual void setup() override;

  virtual void start_of_step() override;
  virtual void post_force() override;
  
  virtual void initial_integrate() override {};

  virtual void pre_final_integrate() override {};
  virtual void final_integrate() override {};
  virtual void post_final_integrate(bool invert_fft = true) override ;
  

  virtual void end_of_step() override {};
private:

  double normalization;
  bool ftphi_flag;
};

}

#endif
