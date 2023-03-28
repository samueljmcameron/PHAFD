#include <algorithm>
#include <iostream>
#include <random>

#include "utility.hpp"


#include "beadrodpmer/initialise.hpp"
#include "beadrodpmer/no_tether.hpp"
#include "beadrodpmer/single_tether.hpp"
#include "beadrodpmer/double_tether.hpp"

#include "atom.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "integrate.hpp"

#include "fixatom_semiflexible.hpp"

using namespace PHAFD_NS;

template <typename T>
FixAtomSemiFlexible<T>::FixAtomSemiFlexible(PHAFD *phafd) : Fix(phafd) {
};


template <typename T>
void FixAtomSemiFlexible<T>::init(const std::vector<std::string> &v_line)
{

  Fix::init(v_line);

  find_group(v_line.at(1));


  // v_line contains a seed in the argument, so need to extract that and initialise
  //  a new RNG with it.
  // need to reconvert to a single string (removing the fix name and the group name)
  std::string line = "";

  std::vector<std::string> new_v_line(v_line);



  int seed = std::stoi(new_v_line.at(2));

  seed = utility::make_unique_seed(seed,world,commbrick->me,commbrick->nprocs);
  
  std::mt19937 gen;
  std::uniform_int_distribution<int> integer_dist;
  
  gen.seed(seed);

  new_v_line.erase(new_v_line.begin(),new_v_line.begin()+3);
  
  for (auto &word : new_v_line)
    line += word + std::string(" ");

  line += "seed ";


  
  
  
  // then add in number of beads from group

  std::string tmpline;


  for (int i = 0; i < start_indices.size(); i++) {
    nbeads.push_back(end_indices[i]-start_indices[i]);
    
    tmpline = std::string("beads ") + std::to_string(nbeads[i])
      + std::string(" ") +  line + std::to_string(integer_dist(gen));
    new_v_line = utility::split_line(tmpline);
    
    
    pmers.push_back(std::make_unique<T>(new_v_line));
  }
  

}

template <typename T>  
void FixAtomSemiFlexible<T>::setup()
{

  int errflag= 0;
  int totalerrflag;


  try {
    check_bonds();
  } catch (const std::runtime_error& excpt) {
    errflag = 1;

  }

  MPI_Allreduce(&errflag,&totalerrflag,1,MPI_INT,MPI_SUM,world);

  if (totalerrflag)
    throw std::runtime_error("Adjacent atom displacements do not match the FixAtomSemiFlexible "
			     "bondlength specified.");
  
  for (int i = 0; i < pmers.size(); i++) {
    pmers[i]->setup(atoms->uxs.middleCols(start_indices[i],nbeads[i]));

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
      
      delta = (atoms->uxs.col(mu)-atoms->uxs.col(mu-1)).norm() - bondlength;
      
      if (abs(delta) > 1e-3) {
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
    pmers[i]->first_step(atoms->uxs.middleCols(start_indices[i],nbeads[i]),
			 atoms->Fs.middleCols(start_indices[i],nbeads[i]),dt);
    for (int bead = start_indices[i]; bead < end_indices[i]; bead++)
      domain->map(atoms->xs.col(bead),atoms->uxs.col(bead),atoms->images[bead]);
  }
  
  
}
template <typename T>
void FixAtomSemiFlexible<T>::final_integrate()
{
  int itermax = 50;

  int iterations;

  int errflag = 0;
  int totalerrflag;
  
  for (int i = 0; i < pmers.size(); i++) {
    iterations =
      pmers[i]->second_step(atoms->uxs.middleCols(start_indices[i],nbeads[i]),
			   atoms->Fs.middleCols(start_indices[i],nbeads[i]),dt,itermax);
    for (int bead = start_indices[i]; bead < end_indices[i]; bead++)
      domain->map(atoms->xs.col(bead),atoms->uxs.col(bead),atoms->images[bead]);

    if (iterations > itermax) {
      std::cerr << "bad timestep from no_tether of polymer on processor "
		<< commbrick->me << " at timestep " << integrate->timestep << std::endl;
      errflag = 1;
    }
    
  }

  MPI_Allreduce(&errflag,&totalerrflag,1,MPI_INT,MPI_SUM,world);

  if (totalerrflag)
    throw std::runtime_error("Failed to maintain polymer bond-length constraints.");


  return;
}


template class FixAtomSemiFlexible<BeadRodPmer::NoTether>;
template class FixAtomSemiFlexible<BeadRodPmer::SingleTether>;
template class FixAtomSemiFlexible<BeadRodPmer::DoubleTether>;
