#include <algorithm>

#include "utility.hpp"


#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "integrate.hpp"
#include "fixgrid_modelb_advect.hpp"
#include "fixgrid_gradient.hpp"
#include "conjugate_noise.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;


FixGridModelBadvect::FixGridModelBadvect(PHAFD *phafd)
  : FixGridModelB(phafd)
{

  for (int dim = 0; dim < 3; dim++)
    flow[dim] = nullptr;
  
  
};


FixGridModelBadvect::~FixGridModelBadvect()
{


}



void FixGridModelBadvect::init(const std::vector<std::string> &v_line,
			       bool add_to_names)
/*
  v_line should have form
  fixname,seed

  followed by the key value pairs

  mobility, value
  temp, value

  in some order.
 */
{

  FixGridModelB::init(v_line,add_to_names);

  std::vector<std::string> new_v_line;
  
  new_v_line.push_back(name+"_phigradient");

  gradphifix = std::make_unique<FixGridGradient>(phafd);
  gradphifix->init(new_v_line,false);

  
  
}


void FixGridModelBadvect::setup()
{

  FixGridModelB::setup();
  gradphifix->setup();

  int Nx = grid->boxgrid[0];
  int Ny = grid->boxgrid[1];
  int Nz = grid->boxgrid[2];

  
  for (int dim = 0; dim < 3; dim++) {
    if (!flow[dim])
      flow[dim] = std::make_unique<fftwArr::array3D<double>>
	(world,"flow[" + std::to_string(dim) + "]" +  name,Nx,Ny,Nz);
  }


  set_flow();
  
}


void FixGridModelBadvect::set_flow()
{

  double y;
  double dy = domain->period[1]/grid->boxgrid[1];

  for (int dim = 0; dim < 3; dim++)
    flow.at(dim)->setZero();
  
  for (int i = 0; i < localNz; i++) {
    
    for (int j = 0; j < localNy; j++) {
      
      y = dy * j - domain->boxlo[1];
      
      for (int k = 0; k < localNx; k++) {
	
	(*flow[0])(i,j,k) = 10.0*(cos(2*M_PI/domain->period[1]*y)+1);
	
      }
    }
  }
}

void FixGridModelBadvect::start_of_step()
{
  FixGridModelB::start_of_step();

  gradphifix->calculate_gradient(grid->ft_phi.get());	
  
}


void FixGridModelBadvect::post_final_integrate(bool invert_fft)
{
  if (invert_fft) {

    fftw_execute(grid->backward_phi);


    for (int i = 0; i < localNz; i++)
      for (int j = 0; j < localNy; j++)
	for (int k = 0; k < localNx; k++) 

	  for (int dim = 0; dim < 3; dim ++)
	    (*grid->phi)(i,j,k) -= (*flow[dim])(i,j,k)
	      *(*gradphifix->gradient[dim])(i,j,k)*dt;
  }
  
  return;
}
