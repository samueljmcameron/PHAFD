#include <algorithm>
#include <iostream>
#include <random>

#include "utility.hpp"



#include "atom.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "integrate.hpp"

#include "fixatom_drag.hpp"

using namespace PHAFD_NS;


FixAtomDrag::FixAtomDrag(PHAFD *phafd) : Fix(phafd) {
};



void FixAtomDrag::init(const std::vector<std::string> &v_line)
{

  Fix::init(v_line);

  find_group(v_line.at(1));
  drag = std::stod(v_line.at(2));


}

  
void FixAtomDrag::setup()
{

}


  
void FixAtomDrag::initial_integrate()
{
  
  
  for (int i = 0; i < start_indices.size(); i++) {  
    for (int j = start_indices[i]; j < end_indices[i]; j++) {
      atoms->xs.col(j) += drag*atoms->Fs.col(j)*dt*0.5;
    }
  }
  
  
}

void FixAtomDrag::final_integrate()
{

  for (int i = 0; i < start_indices.size(); i++) {  
    for (int j = start_indices[i]; j < end_indices[i]; j++) {
      atoms->xs.col(j) += drag*atoms->Fs.col(j)*dt*0.5;
    }
  }
  
  return;
}

