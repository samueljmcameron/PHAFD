#include <algorithm>

#include "utility.hpp"

#include "fixgrid_velocity.hpp"


#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "integrate.hpp"
#include "fixgrid_vdet.hpp"
#include "fixgrid_vtherm.hpp"
#include "conjugate_noise.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;


FixGridVelocity::FixGridVelocity(PHAFD *phafd) : Fix(phafd) {};



void FixGridVelocity::init(const std::vector<std::string> &v_line,
			   bool add_to_names)
/*
  v_line should take form:
  fixname,seedx,seedy,seedz,viscosity,temperature,(immediate)
 */
{

  Fix::init(v_line,add_to_names);

  std::vector<std::string> new_v_line;
  
  new_v_line = v_line;
  if (v_line.size() == 7) {
    if (v_line.at(6) == "immediate")
      immediate_ifft = true;
  } else // vtherm and vdet must have immediate ifft
    new_v_line.push_back("immediate");


  new_v_line.at(0) = name + "_vtherm";
  vtherm = std::make_unique<FixGridVtherm>(phafd);

  vtherm->init(new_v_line,false);

  new_v_line.clear();
  new_v_line.push_back(name + "_vdet");
  new_v_line.push_back(v_line.at(4));
  new_v_line.push_back("immediate");

  vdet = std::make_unique<FixGridVdet>(phafd);
  vdet->init(new_v_line,false);
  
  
}

  
void FixGridVelocity::setup()
{
  vtherm->setup();
  vdet->setup();
}

void FixGridVelocity::reset_dt()
{
  vtherm->reset_dt();
  vdet->reset_dt();
}

void FixGridVelocity::start_of_step()
{
  vdet->start_of_step();
  vtherm->start_of_step();

}


void FixGridVelocity::pre_final_integrate()
{

  
  // compute the two components of velocity
  vdet->pre_final_integrate();
 
  
  // then add them together
  for (int i = 0; i < grid->velocity[0]->Nz(); i++) 
    for (int j = 0; j < grid->velocity[0]->Ny(); j++) 
      for (int k = 0; k < grid->velocity[0]->Nx(); k++) 
	for (int dim = 0; dim < 3; dim++)
	  (*grid->velocity[dim])(i,j,k) =
	    (*grid->vdet[dim])(i,j,k) + (*grid->vtherm[dim])(i,j,k);

}


