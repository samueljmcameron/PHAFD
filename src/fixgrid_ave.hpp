
#ifndef PHAFD_FIXGRID_AVE_HPP
#define PHAFD_FIXGRID_AVE_HPP


#include "fix.hpp"

#include <memory>

namespace PHAFD_NS {

class Compute;

class FixGridAve : public Fix {
public:
  FixGridAve(PHAFD *);

  virtual void init(const std::vector<std::string> &) override;
  
  virtual void setup() override;

  virtual void initial_integrate() override {};
  virtual void post_force() override {};
  virtual void pre_final_integrate() override {};
  virtual void final_integrate() override {};
  virtual void post_final_integrate() override {};
  

  virtual void start_of_step() override;
  virtual void end_of_step() override;

private:

  int every, repeat, freq;

  int step_counter,next_freq;

  std::vector<int> savesteps;

  bool collect_this_step;
  
  Compute *compute;
  Fix *fix;

  
};

}

#endif
