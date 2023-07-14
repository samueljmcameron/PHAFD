
#ifndef PHAFD_FIXGRID_MODELB_HPP
#define PHAFD_FIXGRID_MODELB_HPP


#include "fix.hpp"

#include <memory>

namespace PHAFD_NS {
class FixGridGradPhi;
class ConjugateVolFrac;

class FixGridModelB : public Fix {
public:
  FixGridModelB(PHAFD *);

  virtual void init(const std::vector<std::string> &) override;
  
  virtual void setup() override;

  virtual void start_of_step() override;
  
  virtual void initial_integrate() override {};
  virtual void post_force() override {};
  virtual void pre_final_integrate() override;
  virtual void final_integrate() override;
  virtual void post_final_integrate() override;
  
  virtual void reset_dt() override;

  virtual void end_of_step() override {};
private:

  std::unique_ptr<ConjugateVolFrac> conjugate;

  bool didnotintegrate;

  std::unique_ptr<FixGridGradPhi> fixgridgradphi;
  
};

}

#endif
