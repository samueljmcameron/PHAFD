#ifndef PHAFD_FIXGRID_HPP
#define PHAFD_FIXGRID_HPP

#include <vector>
#include <string>
#include "pointers.hpp"

namespace PHAFD_NS {

class FixGrid : protected Pointers {
public:
  FixGrid(PHAFD *);

  virtual void init(const std::vector<std::string> &) ;
  
  virtual void setup() = 0;
  virtual void initial_integrate() = 0;
  virtual void final_integrate() = 0;
  


  virtual void reset_dt(double);

  inline static std::vector<std::string> NAMES;

  std::string name;
protected:
  // group properties, all groups must be in chunks

  double dt;

private:

};


}

#endif
