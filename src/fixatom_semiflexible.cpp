#include <iostream>


#include "input.hpp"

#include "beadrodpmer/initialise.hpp"
#include "beadrodpmer/no_tether.hpp"
#include "beadrodpmer/single_tether.hpp"
#include "beadrodpmer/double_tether.hpp"


#include "fixatom_semiflexible.hpp"

using namespace PHAFD;

template <typename T>
FixAtomSemiFlexible<T>::FixAtomSemiFlexible(PHAFD *phafd) : Pointers(phafd) {};



void FixAtomSemiFlexible<T>::init(const std::vector<std::string> &v_line)
{



  FixAtom::init(v_line);

  // need to reconvert to a single string (removing the fix name and the group name)
  std::string line;
  for (std::vector<std::string>::iterator it = v_line.begin()+2; it !=v_line.end(); it++)
    line += *it + std::string(" ");



  // then add in number of beads from group

  std::string tmpline;
  std::vector<std::string> new_v_line;

  for (int i = 0; i < start_indices.size(); i++) {
    nbeads.push_back(end_indices[i]-start_indices[i]);
    
    tmpline = std::string("beads ") + std::to_string(nbeads[i])
      + std::string(" ") +  line;
    new_v_line = input::split_line(tmpline);
    replace_with_new_seed(new_v_line);

    pmers.push_back(std::make_unique<T>(new_v_line));
  }
  

}



template <typename T>  
void FixAtomSemiFlexible<T>::setup()
{

  check_bonds();
  
  for (int i = 0; i < pmers.size(); i++) {
    pmers[i]->setup(atoms.uxs.middleCols(start_indices[i],nbeads[i]));

  }

}



template <typename T>
void FixAtomSemiFlexible<T>::check_bonds()
{
  
  
  double delta;
  double bondlength;
  for (int i = 0; i < pmers.size(); i++) {
    
    bondlength = pmers[i]->get_bondlength();
    
    for (int mu = start_indices[i]+1; mu < end_indices[i]; mu++) {
      
      delta = (atoms.uxs.col(mu)-atoms.uxs.col(mu-1)).norm() - bondlength;
      
      if (abs(delta) > 1e-3) {
	std::cout << "on id " << domain.me << " delta[" << mu << "] = " << delta << std::endl;
	
	throw std::runtime_error("Adjacent atom displacements do not match the FixAtomSemiFlexible "
				 "bondlength specified.");

      }
      
    }
    
    
  } 
    
}



template <typename T>  
void FixAtomSemiFlexible<T>::initial_integrate()
{
  
  
  
  for (int i = 0; i < pmers.size(); i++) {
    pmers[i]->first_step(atoms.uxs.middleCols(start_indices[i],nbeads[i]),
			 atoms.Fs.middleCols(start_indices[i],nbeads[i]),dt);
    for (int bead = start_indices[i]; bead < end_indices[i]; bead++)
      domain.map(atoms.xs.col(bead),atoms.uxs.col(bead),atoms.images[bead]);
  }
  
  
}
template <typename T>
void FixAtomSemiFlexible<T>::final_integrate()
{
  int itermax = 50;

  int iterations;
  
  for (int i = 0; i < pmers.size(); i++) {
    iterations =
      pmers[i]->second_step(atoms.uxs.middleCols(start_indices[i],nbeads[i]),
			   atoms.Fs.middleCols(start_indices[i],nbeads[i]),dt,itermax);
    for (int bead = start_indices[i]; bead < end_indices[i]; bead++)
      domain.map(atoms.xs.col(bead),atoms.uxs.col(bead),atoms.images[bead]);

    if (iterations > itermax)
      std::cerr << "bad timestep from no_tether." << std::endl;
  }


  return;
}


/* given a vector like {"alpha", "beta", "seed", "8492109"}, generate
   a new seed from "8492109" which is unique to the current MPI process.*/
void FixAtomSemiFlexible::replace_with_new_seed(std::string &v_line)
{


  auto pos = std::find(v_line.begin(),v_line.end(),"seed");

  if (pos == v_line.end())
    throw std::runtime_error("No seed variable present in fix/semiflexible.");

  pos ++;

  if (pos == v_line.size())
    throw std::runtime_error("Invalid seed argument in fix/semiflexible.");


  int seed = std::stod(v_line[pos]);

  psPDE::RandomPll rpll(world,commbrick->me,seed,commbrick->nprocs);

  seed = rpll.get_processor_seed();

  v_line[pos] = std::to_string(seed);

}



template class FixAtomSemiFlexible<BeadRodPmer::NoTether>;
template class FixAtomSemiFlexible<BeadRodPmer::SingleTether>;
template class FixAtomSemiFlexible<BeadRodPmer::DoubleTether>;
