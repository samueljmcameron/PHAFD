#ifndef PHAFD_GROUP_HPP
#define PHAFD_GROUP_HPP

#include <vector>
#include <string>

#include "pointers.hpp"

namespace PHAFD_NS {

class Group : protected Pointers {
public:
  Group(PHAFD *);
  
  void create_all();
  void create_group(const std::vector<std::string> &);

  // iterate from start_indices[i] to end_indices[i]
  //  - length of these two vectors will be 1 if atom group,
  //  and >= 1 if molecule group

  
  std::vector<int> start_indices; // first index in chunk
  std::vector<int> end_indices; // last index + 1 in chunk

  std::string name,style;
  
private:
  // group properties, all groups must be in chunks

  void group_atoms(int,int);
  void group_molecule(int);

  
};


}
#endif
