#include "group.hpp"

#include "input.hpp"
#include <stdexcept>
#include <set>

using namespace PHAFD;


Group::Group(Atom &atoms)
{
  name = "all";
  style = "atom";

  start_indices.push_back(0);
  end_indices.push_back(atoms.nowned);


}

Group::Group(std::string line, Atom &atoms)
{


  std::set<int> group_set;

  std::vector<std::string> v_line = input::split_line(line);

  name = v_line[0];

  if (name == "all")
    throw std::invalid_argument("Cannot call group reserved word 'all'.");
  
  style = v_line[1];
  
  if (v_line[2] == "<>") {

    if (v_line.size() != 5)
      throw std::invalid_argument("need start and end point for group using '<>' argument.");

    int nstart = std::stoi(v_line[3]);
    int nend = std::stoi(v_line[4]);

    if (nend < nstart)
      throw std::invalid_argument("end point for group larger than start point when using '<>'.");
    if (nstart < 0)
      throw std::invalid_argument("start point of group must be larger than zero.");    
    for (int i = nstart; i <= nend; i++)
      group_set.insert(i);
    
    
  } else {

    if (v_line.size() < 3)
      throw std::invalid_argument("No atoms/molecules specified for group.");


    int tag;
    for (std::string::size_type i = 2; i < v_line.size(); i++) {
      tag = std::stoi(v_line[i]);
      if (tag < 0)
	throw std::invalid_argument("Group atom/molecule IDs must be greater than zero.");
      group_set.insert(tag);
    }
    
  }


  
  if (style == "molecule") {

    for (auto imol : group_set)
      group_molecule(imol,atoms);
    

  } else if (style == "atom") {

    if (v_line[2] != "<>")
      throw std::invalid_argument("Group atom requires '<>' style parameters.");

    group_atoms(*group_set.begin(),*group_set.rbegin(),atoms);


  } else {
    throw std::invalid_argument("Group must be of style molecule or atom.");
  }

}


void Group::group_atoms(int tagstart, int tagend, const Atom &atoms)
{


  if (tagstart < 0 || tagend < tagstart)
    throw std::invalid_argument("Invalid IDs in group atoms.");

  

  int startatom = 0;

  // loop through atoms to see if the first element of the set is in the atoms
  while (startatom < atoms.nowned && atoms.tags[startatom] < tagstart) {
    startatom ++;

  }



  int tag = tagstart;


  // if the very first atom ID is larger than (or equal to) the first ID of the group

  while (tag <= tagend && atoms.tags[startatom] != tag) tag ++;

    
  if (tag <= tagend)     // if this is true, then atom is part of the group
    start_indices.push_back(startatom);
  else // otherwise, the atom IDs are all bigger than all members of the group
    return;


  
  int iatom = startatom;

  while (iatom < atoms.nowned && tag <= tagend) {
    if (atoms.tags[iatom] != tag) 
      throw std::invalid_argument("atoms of the same group must be adjacently owned.");

    iatom ++;
    tag ++;
    
  }
  end_indices.push_back(iatom);
    
  
  return;

}

void Group::group_molecule(int imol, const Atom &atoms)
{

  int iatom = 0;
  int count = 0;
  int lastindex;


  // loop over all atoms, see whether atoms are stored correctly and if molecule is on processor.
  while (iatom < atoms.nowned) {
    
    if (atoms.moltags[iatom] == imol ) {
      
      if (count == 0)
	start_indices.push_back(iatom);
      else if (lastindex != iatom-1) 
	throw std::runtime_error("atoms of the same molecule must be adjacently owned.");

      count += 1;
      lastindex = iatom;
    }
    
    iatom += 1;
  }


  if (count != 0)
    end_indices.push_back(lastindex+1);

  
  
}
