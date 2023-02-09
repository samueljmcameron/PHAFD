#ifndef PHAFD_INTEGRATE_HPP
#define PHAFD_INTEGRATE_HPP

#include "pointers.hpp"

namespace PHAFD_NS {

class Integrate : protected Pointers {
public:
  Integrate(PHAFD *);

  ~Integrate();

  void setup();
  void run();
  void run_until_touching(double);

  int nsteps; // total number of steps
  int firststep; // first step
  int timestep; // current step
  double dt; // time interval between steps
  
};
  
}
#endif
