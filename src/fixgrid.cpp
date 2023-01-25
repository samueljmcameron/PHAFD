
#include "group.hpp"

#include "fixgrid.hpp"

using namespace PHAFD_NS;


FixGrid::FixGrid(PHAFD *phafd) : Pointers(phafd) {};


void FixGrid::init(const std::vector<std::string> &v_line) 
{

  if (v_line.size() < 2)
    throw std::runtime_error("incorrect args in fixgrid.");
    
  name = v_line.at(0);


  for (auto &fixname : NAMES)
    if (fixname == name)
      throw std::runtime_error("Error: Duplicate of fix " + name + std::string("."));

  NAMES.push_back(name);


}



void FixGrid::reset_dt(double timestep)
{

  dt = timestep;
}

