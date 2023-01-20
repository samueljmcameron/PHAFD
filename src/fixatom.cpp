#include "fixatom.hpp"


using namespace PHAFD;


FixAtom::FixAtom(const Group &group, std::string line,Atom &atoms,Domain &domain)
  : start_indices(group.start_indices),end_indices(group.end_indices),atoms(atoms),
    domain(domain)
{

}


void FixAtom::reset_dt(double timestep)
{

  dt = timestep;
}
