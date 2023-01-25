
#include "group.hpp"

#include "fixatom.hpp"

using namespace PHAFD_NS;


FixAtom::FixAtom(PHAFD *phafd) : Pointers(phafd) {};


void FixAtom::init(const std::vector<std::string> &v_line) 
{

  if (v_line.size() < 3)
    throw std::runtime_error("incorrect args in fix/semiflexible.");
    
  name = v_line.at(0);


  for (auto &fixname : NAMES)
    if (fixname == name)
      throw std::runtime_error("Error: Duplicate of fix " + name + std::string("."));

  NAMES.push_back(name);

  // sets start_indices and end_indices
  find_group(v_line.at(1));

}


void FixAtom::find_group(const std::string &gname)
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


void FixAtom::reset_dt(double timestep)
{

  dt = timestep;
}

