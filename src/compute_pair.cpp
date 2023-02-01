#include "compute_pair.hpp"

#include <cmath>

#include "pair.hpp"
#include "atom.hpp"
#include "comm_brick.hpp"

using namespace PHAFD_NS;


ComputePair::ComputePair(PHAFD *phafd) : Compute(phafd) {
  per_atom = true;
  numberofcomponents = 3;
}


void ComputePair::init(const std::vector<std::string> &v_line) {

  Compute::init(v_line);


  std::string pairname = v_line.at(1);

  pair = nullptr;  

  for (auto &ps: pairs) {
    if (pairname == ps->name) {
      pair = ps.get();
      break;
    }
  }

  if (pair == nullptr)
    throw std::invalid_argument("Cannot compute pair style " + pairname
				+ ", it hasn't been created!");


  array.resize(atoms->nowned*numberofcomponents);
  
}

void ComputePair::end_of_step()
{

  // store total atom forces temporarily in the array, and zero the atom forces
  int index = 0;
  for (int iatom = 0; iatom < atoms->nowned; iatom++) {
    array[index++] = atoms->Fs(0,iatom);
    array[index++] = atoms->Fs(1,iatom);
    array[index++] = atoms->Fs(2,iatom);

  }

  atoms->Fs.setZero();
  // now compute the forces due exclusively to this pair style, and put them
  // in atoms->Fs
  pair->compute();
  commbrick->reverse_comm();

  // now swap the values of array and Fs, so that Fs ends up having the cumulative
  // forces on the atoms (as it did before calling this compute), and now the
  // array has the forces on the atoms due to this pair style exclusively
  index = 0;
  double tmp;
  for (int iatom = 0; iatom < atoms->nowned; iatom++) {

    // store cumulative force in a temporary
    tmp = array[index];
    // then put pair force into array
    array[index++] = atoms->Fs(0,iatom);
    // finally put cumulative force back where it belongs
    atoms->Fs(0,iatom) = tmp;

    tmp = array[index];
    array[index++] = atoms->Fs(1,iatom);
    atoms->Fs(1,iatom) = tmp;


    tmp = array[index];
    array[index++] = atoms->Fs(2,iatom);
    atoms->Fs(2,iatom) = tmp;
    
  }
}
