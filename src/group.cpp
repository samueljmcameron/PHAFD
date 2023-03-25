#include "group.hpp"
#include <stdexcept>
#include <set>

#include "atom.hpp"

using namespace PHAFD_NS;


Group::Group(PHAFD *phafd) : Pointers(phafd)
{};

void Group::create_all()
{
  name = "all";
  style = "atom";


  for (auto &gname : NAMES)
    if (gname == name)
      throw std::runtime_error("Error: Duplicate of group " + name + std::string("."));

  NAMES.push_back(name);
  
  start_indices.push_back(0);
  end_indices.push_back(atoms->nowned);

}


void Group::create_group(const  std::vector<std::string> &v_line)
{


  std::set<int> group_set;


  name = v_line.at(0);

  for (auto &gname : NAMES)
    if (gname == name)
      throw std::runtime_error("Error: Duplicate of group " + name + std::string("."));

  NAMES.push_back(name);

  

  if (name == "all")
    throw std::invalid_argument("Cannot call group reserved word 'all'.");
  
  style = v_line.at(1);
  
  if (v_line.at(2) == "<>") {

    if (v_line.size() != 5)
      throw std::invalid_argument("need start and end point for group using '<>' argument.");

    int nstart = std::stoi(v_line.at(3));
    int nend = std::stoi(v_line.at(4));

    if (nend < nstart)
      throw std::invalid_argument("end point for group larger than start point when using '<>'.");
    if (nstart < 0)
      throw std::invalid_argument("start point of group must be larger than zero.");    
    for (int i = nstart; i <= nend; i++)
      group_set.insert(i);
    
    
  } else {

    if (v_line.size() < 3)
      throw std::invalid_argument("No atoms/molecules specified for group.");


    std::size_t pos;
    int tag;
    for (std::string::size_type i = 2; i < v_line.size(); i++) {
      tag = std::stoi(v_line.at(i),&pos);
      if (pos != v_line.at(i).length()) 
	throw std::invalid_argument("Group atom/molecule ID is not an integer.");
      if (tag < 0)
	throw std::invalid_argument("Group atom/molecule IDs must be greater than zero.");
      group_set.insert(tag);
    }
    
  }


  
  if (style == "molecule") {

    for (auto imol : group_set)
      group_molecule(imol);

  } else if (style == "atom") {

    if (v_line.at(2) != "<>")
      throw std::invalid_argument("Group atom requires '<>' style parameters.");

    group_atoms(*group_set.begin(),*group_set.rbegin());


  } else {
    throw std::invalid_argument("Group must be of style molecule or atom.");
  }

}


void Group::group_atoms(int tagstart, int tagend)
{


  if (tagstart < 0 || tagend < tagstart)
    throw std::invalid_argument("Invalid IDs in group atoms.");

  

  int startatom = 0;

  // loop through atoms to see if the first element of the set is in the atoms
  while (startatom < atoms->nowned && atoms->tags[startatom] < tagstart) {
    startatom ++;

  }



  int tag = tagstart;
  int foundatoms = 0;

  // if the very first atom ID is larger than (or equal to) the first ID of the group

  while (tag <= tagend && atoms->tags[startatom] != tag) tag ++;

    
  if (tag <= tagend) {     // if this is true, then atom is part of the group
    foundatoms = 1;
    start_indices.push_back(startatom);

    

  
    int iatom = startatom;
    
    while (iatom < atoms->nowned && tag <= tagend) {
      if (atoms->tags[iatom] != tag) 
	throw std::invalid_argument("atoms of the same group must be adjacently owned.");
      
      iatom ++;
      tag ++;
      
    }
    end_indices.push_back(iatom);
  }
    
  int all_foundatoms;

  MPI_Allreduce(&foundatoms,&all_foundatoms,1,MPI_INT,MPI_SUM,world);

  if (all_foundatoms == 0)
    throw std::runtime_error("atoms " + std::to_string(tagstart) + std::string(" to ")
			     + std::to_string(tagend) 
			     + std::string(" do not exist on any processor." ));

  return;

}

void Group::group_molecule(int imol)
{

  int iatom = 0;
  int count = 0;
  int lastindex;


  // loop over all atoms, see whether atoms are stored correctly and if molecule is on processor.
  while (iatom < atoms->nowned) {
    
    if (atoms->moltags[iatom] == imol ) {
      
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

  int totalcount;

  MPI_Allreduce(&count,&totalcount,1,MPI_INT,MPI_SUM,world);

  if (totalcount == 0)
    throw std::runtime_error("Group error: Molecule " + std::to_string(imol)
			     + std::string(" does not exist on any processor." ));
  
}
