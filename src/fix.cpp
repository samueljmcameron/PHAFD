
#include "group.hpp"
#include "fix.hpp"

#include "integrate.hpp"

using namespace PHAFD_NS;


Fix::Fix(PHAFD *phafd) : Pointers(phafd) {
  

  per_grid = per_ftgrid = per_atom = scalar = vector = false;
  averaging = false;
  numberofcomponents = 1;
  num_less_zero = num_great_zero = 0;
};


void Fix::init(const std::vector<std::string> &v_line) 
{
  this_step = false;
  if (v_line.size() < 1)
    throw std::runtime_error("incorrect args in fixgrid.");
    
  name = v_line.at(0);


  for (auto &fixname : NAMES)
    if (fixname == name)
      throw std::runtime_error("Error: Duplicate of fix " + name + std::string("."));

  NAMES.push_back(name);


}



void Fix::find_group(const std::string &gname)
{
  int group_index = -1;
  
  int counter = 0;
  for (auto &groupname : Group::NAMES) {
    if (groupname == gname)
      group_index = counter;

    counter ++;
  }

  if (group_index == -1)
    throw std::runtime_error("Group " + gname + std::string(" does not exist for fix ")
			     + name);


  start_indices = groups.at(group_index)->start_indices;
  end_indices = groups.at(group_index)->end_indices;

}




void Fix::reset_dt()
{

  dt = integrate->dt;
}

void Fix::start_of_step()
{
  this_step = false;
}
