
#ifndef PHAFD_FIXGRID_FLORYHUGGINS_HPP
#define PHAFD_FIXGRID_FLORYHUGGINS_HPP


#include "fix.hpp"

#include <memory>

namespace psPDE {
  class FixGridFloryHuggins;
}


namespace PHAFD_NS {


class FixGridFloryHuggins : public Fix {
public:
  FixGridFloryHuggins(PHAFD *);

  virtual void init(const std::vector<std::string> &) override;
  
  virtual void setup() override;

  virtual void initial_integrate() override {};
  virtual void post_force() override;
  virtual void pre_final_integrate() override {};
  virtual void final_integrate() override {};
  virtual void post_final_integrate() override {};

  virtual void end_of_step() override {};
  
private:

  std::unique_ptr<psPDE::FixGridFloryHuggins> ps_flory;
  
};

}

#endif
