
#ifndef PHAFD_FIXGRID_CONJUGATE_HPP
#define PHAFD_FIXGRID_CONJUGATE_HPP


#include "fixgrid.hpp"

#include <memory>

namespace PHAFD_NS {

template <typename T>
class FixGridConjugate : public FixGrid {
public:
  FixGridConjugate(PHAFD *);

  virtual void init(const std::vector<std::string> &) override;
  
  virtual void setup() override;

  virtual void initial_integrate() override;

  virtual void pre_final_integrate() override {};
  virtual void final_integrate() override;

  virtual void reset_dt(double) override;
private:

  std::unique_ptr<T> conjugate;

  
};

}

#endif
