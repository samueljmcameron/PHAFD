#ifndef PHAFD_INTEGRATE_HPP
#define PHAFD_INTEGRATE_HPP
#include <cstdint>
#include "pointers.hpp"

namespace PHAFD_NS {

class Integrate : protected Pointers {
public:
  Integrate(PHAFD *);

  ~Integrate();

  void setup();
  void run();
  void run_until_touching(double);

  int64_t nsteps; // total number of steps
  int64_t firststep; // first step
  int64_t timestep; // current step
  double dt; // time interval between steps
private:
  void finalise();
};
  
}
#endif
