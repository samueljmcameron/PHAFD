
#ifndef PHAFD_FIXGRID_FLORYHUGGINS_HPP
#define PHAFD_FIXGRID_FLORYHUGGINS_HPP


#include "fix.hpp"

#include <memory>



namespace PHAFD_NS {


class FixGridFloryHuggins : public Fix {
public:
  FixGridFloryHuggins(PHAFD *);

  virtual void init(const std::vector<std::string> &,bool) override;
  
  virtual void setup() override;

  virtual void initial_integrate() override {};
  virtual void post_force() override;
  virtual void pre_final_integrate() override {};
  virtual void final_integrate() override {};
  virtual void post_final_integrate(bool /* flag */) override {};

  virtual void end_of_step() override {};
  
private:

  double chemical_potential(double);

  double temp,volFH,chi,kappa;
};

}

#endif
