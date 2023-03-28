
#ifndef PHAFD_FIXATOM_DRAG_HPP
#define PHAFD_FIXATOM_DRAG_HPP


#include "fix.hpp"

#include <memory>

namespace PHAFD_NS {

class FixAtomDrag : public Fix {
public:
  FixAtomDrag(PHAFD *);

  virtual void init(const std::vector<std::string> &) override;
  
  virtual void setup() override;
  virtual void initial_integrate() override;
  virtual void post_force() override {};
  virtual void pre_final_integrate() override {};
  virtual void final_integrate() override;
  virtual void post_final_integrate() override {};

  virtual void end_of_step() override {};

private:

  double drag;

  
};

}

#endif
