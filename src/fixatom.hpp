#ifndef PHAFD_FIXATOM_HPP
#define PHAFD_FIXATOM_HPP

#include <vector>
#include <string>
#include "atom.hpp"
#include "group.hpp"
#include "domain.hpp"

namespace PHAFD {

class FixAtom {
public:
  FixAtom(const Group &,std::string,Atom &,Domain &);

  virtual void setup() = 0;
  virtual void initial_integrate() = 0;
  virtual void final_integrate() = 0;

  virtual void reset_dt(double);

protected:
  // group properties, all groups must be in chunks

  const std::vector<int> start_indices,end_indices;

  Atom &atoms;
  Domain &domain;
  double dt;
};


}

#endif
