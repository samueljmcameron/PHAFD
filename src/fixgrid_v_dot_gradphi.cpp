#include <algorithm>

#include "utility.hpp"

#include "fixgrid_v_dot_gradphi.hpp"


#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "integrate.hpp"
#include "fixgrid_velocity.hpp"
#include "fixgrid_vdet.hpp"
#include "fixgrid_vtherm.hpp"
#include "conjugate_noise.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;


FixGridVdotGradPhi::FixGridVdotGradPhi(PHAFD *phafd) : Fix(phafd) {};



void FixGridVdotGradPhi::init(const std::vector<std::string> &v_line,
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


  new_v_line.at(0) = name + "_velocity";
  velocity = std::make_unique<FixGridVelocity>(phafd);
  velocity->init(new_v_line,false);
  
}

  
void FixGridVdotGradPhi::setup()
{
  velocity->setup();
  normalization = (grid->boxgrid[0]*grid->boxgrid[1]*grid->boxgrid[2]);


}

void FixGridVdotGradPhi::reset_dt()
{
  velocity->reset_dt();
}

void FixGridVdotGradPhi::start_of_step()
{
  velocity->start_of_step();

}


void FixGridVdotGradPhi::pre_final_integrate()
{

  
  // compute the two components of velocity
  velocity->pre_final_integrate();
 
  // then add them together
  for (int i = 0; i < grid->v_dot_gradphi->Nz(); i++) 
    for (int j = 0; j < grid->v_dot_gradphi->Ny(); j++) 
      for (int k = 0; k < grid->v_dot_gradphi->Nx(); k++) {
	(*grid->v_dot_gradphi)(i,j,k) = 0;
	for (int dim = 0; dim < 3; dim++)
	  (*grid->v_dot_gradphi)(i,j,k) +=
	    (*grid->velocity[dim])(i,j,k)*(*grid->gradphi[dim])(i,j,k);

      }


  if (!immediate_ifft) 
    fftw_execute(grid->forward_v_dot_gradphi);

}


void FixGridVdotGradPhi::post_final_integrate(bool invert_fft)
{
  if (immediate_ifft)
    return;
  if (invert_fft) {
    fftw_execute(grid->backward_v_dot_gradphi);
    (*grid->v_dot_gradphi) /= normalization;
  }

}


