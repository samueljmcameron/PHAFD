#include "input.hpp"


#include "utility.hpp"


#include "domain.hpp"
#include "grid.hpp"
#include "group.hpp"
#include "comm_brick.hpp"
#include "read_atoms.hpp"
#include "atom.hpp"




#include <iostream>
#include <string>
#include <set>

using namespace PHAFD_NS;

Input::Input(PHAFD *phafd) : Pointers(phafd) {}

Input::~Input() = default;

void Input::read()
{
  std::vector<std::string> v_line;
  std::string line;

  line = "boxdims 100.0 100.0 100.0 boxorigin -50.0 -50.0 -50.0";
  v_line = utility::split_line(line);
  domain->set_box(v_line);


  line = "64 64 64 concentration";
  v_line = utility::split_line(line);

  grid->create(v_line);
  domain->set_subbox();

  line = "constant concentration 0.2 0.2 891409";
  v_line = utility::split_line(line);
  grid->populate(v_line);


  std::string atom_data_fname = "read_p%.data";

  utility::replacePercentages(atom_data_fname,commbrick->me);


  

  ReadAtoms read_atoms(phafd);

  read_atoms.read_file(atom_data_fname);

  atoms->check_tags_and_types();

  
  groups.push_back(std::make_unique<Group>(phafd));
  groups[0]->create_all();



  groups.push_back(std::make_unique<Group>(phafd));
  line = "nt_pmer0 molecule 0";
  v_line = utility::split_line(line);    
  groups[1]->create_group(v_line);


  groups.push_back(std::make_unique<Group>(phafd));
  line = "nt_pmer1 molecule 1";
  v_line = utility::split_line(line);    
  groups[2]->create_group(v_line);

  groups.push_back(std::make_unique<Group>(phafd));
  line = "nt_pmer2 molecule 2";
  v_line = utility::split_line(line);    
  groups[3]->create_group(v_line);

  groups.push_back(std::make_unique<Group>(phafd));
  line = "nt_pmer3 molecule 3";
  v_line = utility::split_line(line);    
  groups[4]->create_group(v_line);

  std::set<std::string> gnames;
  for (auto &group : groups) {
    gnames.insert(group->name);
  }

  if (gnames.size() != groups.size())
    throw std::runtime_error("Cannot have two groups with same name.");
    

  line = "bondlength 2.0 zeta_para 0.5 zeta_perp 1.0 bending 206.0 temp 4.114 seed "
      + std::to_string(seed0);
  fixes.push_back(std::make_unique<Fix>(phafd));
  fixes[0].init()
  

  
  if (commbrick->me == 0)
    std::cout << groups[0]->name << std::endl;
  
}
