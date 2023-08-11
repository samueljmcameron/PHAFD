#include "compute.hpp"

using namespace PHAFD_NS;


Compute::Compute(PHAFD *phafd) : Pointers(phafd) {
  per_grid = per_ftgrid = per_atom = scalar = vector =  false;
  clusterscomputed = false;
  numberofcomponents = 1;
}


void Compute::init(const std::vector<std::string> &v_line) {

  this_step = false;
  
  name = v_line.at(0);

  for (auto &computename : NAMES)
    if (computename == name)
      throw std::runtime_error("Error: Duplicate of compute " + name + std::string("."));

  NAMES.push_back(name);
  

}

void Compute::start_of_step()
{
  this_step = false;
}
