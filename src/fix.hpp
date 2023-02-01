#ifndef PHAFD_FIX_HPP
#define PHAFD_FIX_HPP

#include <vector>
#include <string>
#include <set>
#include "pointers.hpp"

namespace PHAFD_NS {

class Fix : protected Pointers {
public:
  Fix(PHAFD *);

  virtual void init(const std::vector<std::string> &) ;


  // must call setup AFTER dump setup
  virtual void setup() = 0;
  virtual void initial_integrate() = 0;
  virtual void post_force() = 0;
  virtual void pre_final_integrate() = 0;
  virtual void final_integrate() = 0;
  virtual void post_final_integrate() = 0;

  virtual void reset_dt();
  virtual void start_of_step();
  virtual void end_of_step() = 0;
  
  
  inline static std::vector<std::string> NAMES;

  std::string name;
  bool per_grid;
  bool per_ftgrid;
  bool per_atom;
  bool global;

  bool this_step;

  std::set<std::string> dump_callers;
  std::vector<double> array;

  int Nx,Ny,Nz;
  int numberofcomponents; // same as compute
  
protected:
  // group properties, all groups must be in chunks, only relevant for atom fixes at the moment
  std::vector<int> start_indices,end_indices;
  void find_group(const std::string &)  ;
  double dt;

};


}

#endif
