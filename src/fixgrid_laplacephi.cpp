#include <algorithm>
#include <iostream>
#include "utility.hpp"

#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "integrate.hpp"
#include "fixgrid_laplacephi.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;


FixGridLaplacePhi::FixGridLaplacePhi(PHAFD *phafd) : Fix(phafd) {
  once = false;
};



void FixGridLaplacePhi::init(const std::vector<std::string> &v_line)
{

  Fix::init(v_line);

  if (v_line.size() == 2) {
    if (v_line.at(1) == "once") once = true;
    else if (v_line.at(1) == "every") once = false;
  } else {
    throw std::runtime_error("need to specify when to compute gradient in fix/laplacephi");
  }
  
}

  
void FixGridLaplacePhi::setup()
{

}


void FixGridLaplacePhi::start_of_step()
{


  if ((once == true && integrate->timestep == integrate->firststep+1) || !once) {
    
    // calculate the laplacian of phi.
    fftw_execute(grid->forward_phi);
    
    const int local0start = grid->ft_phi->get_local0start();
    const int globalNy = grid->ft_boxgrid[1];
    const int globalNz = grid->ft_boxgrid[2];
    
    const double dqy = domain->dqy();
    const double dqz = domain->dqz();
    const double dqx = domain->dqx();


    double l,m,n;
    double neg_q2;
    
    for (int i = 0; i < grid->ft_phi->Nz(); i++) {
      
      if (i + local0start > globalNz/2) 
	l = -globalNz + i + local0start;
      else
	l = i + local0start;
      
      for (int j = 0; j < grid->ft_phi->Ny(); j++) {
	
	if (j > globalNy/2)
	  m = -globalNy + j;
	else
	  m = j;
	
	for (int k = 0; k < grid->ft_phi->Nx(); k++) {
	  
	  n = k;

	  neg_q2 = -(dqx*n*dqx*n + dqy*m*dqy*m + dqz*l*dqz*l);

	  
	  (*grid->ft_laplacephi)(i,j,k) = (*grid->ft_phi)(i,j,k)*neg_q2;
;
	  
	}
      }
    }
    
    fftw_execute(grid->backward_laplacephi);
    
    double factor = grid->boxgrid[0]*grid->boxgrid[1]*grid->boxgrid[2];
    *grid->laplacephi /= factor;
    
  }


}


