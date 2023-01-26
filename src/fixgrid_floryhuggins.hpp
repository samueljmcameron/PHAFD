
#ifndef PHAFD_FIXGRID_FLORYHUGGINS_HPP
#define PHAFD_FIXGRID_FLORYHUGGINS_HPP


#include "fixgrid.hpp"

#include <memory>

namespace psPDE {
  class FixGridFloryHuggins;
}


namespace PHAFD_NS {


class FixGridFloryHuggins : public FixGrid {
public:
  FixGridFloryHuggins(PHAFD *);

  virtual void init(const std::vector<std::string> &) override;
  
  virtual void setup() override;

  virtual void initial_integrate() override {};
  virtual void pre_final_integrate() override;
  virtual void final_integrate() override {};

private:

  std::unique_ptr<psPDE::FixGridFloryHuggins> ps_flory;
  
};

}

#endif
