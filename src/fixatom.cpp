#include "fixatom.hpp"


using namespace PHAFD_NS;


FixAtom::FixAtom(PHAFD *phafd) : Pointers(phafd) {};


void FixAtom::init(const std::vector<std::string> &v_line) 
{

  if (v_line.size() < 3)
    throw std::runtime_error("incorrect args in fix/semiflexible.");
    
  name = v_line[0];

  // sets start_indices and end_indices
  find_group(v_line[1]);

}


void FixAtom::find_group(const std::string &gname)
{
  group_index = -1;
  
  int counter;
  for (auto & group : groups) {
    if (group->name == gname)
      group_index = counter;

    counter ++;
  }

  if (group_index == -1)
    throw std::runtime_error("Group " + group_name + std::string(" does not exist for fix ")
			     + name);


  start_indices = groups[group_index]->start_indices;
  end_indices = groups[group_index]->end_indices;

}


void FixAtom::reset_dt(double timestep)
{

  dt = timestep;
}
