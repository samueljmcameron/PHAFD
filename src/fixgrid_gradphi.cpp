#include <algorithm>

#include "utility.hpp"

#include "ps_pde/grid.hpp"
#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "integrate.hpp"
#include "fixgrid_gradphi.hpp"

using namespace PHAFD_NS;


FixGridGradPhi::FixGridGradPhi(PHAFD *phafd) : Fix(phafd) {
  once = false;
};



void FixGridGradPhi::init(const std::vector<std::string> &v_line)
{

  Fix::init(v_line);

  if (v_line.size() == 2) {
    if (v_line.at(1) == "once") once = true;
    else if (v_line.at(1) == "every") once = false;
  } else {
    throw std::runtime_error("need to specify when to compute gradient in fix/gradphi");
  }
  
}

  
void FixGridGradPhi::setup()
{

}


void FixGridGradPhi::start_of_step()
{

  if ((once == true && integrate->timestep == integrate->firststep) || !once) {
    
    // calculate the gradient of phi.
    fftw_execute(grid->ps_grid->forward_phi);
    
    const int local0start = grid->ft_phi->get_local0start();
    const int globalNy = grid->ft_boxgrid[1];
    const int globalNz = grid->ft_boxgrid[2];
    
    
    std::complex<double> idqx(0,domain->dqx());
    std::complex<double> idqy(0,domain->dqy());
    std::complex<double> idqz(0,domain->dqz());
    
    double l,m,n;
    
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
	  
	  (*grid->ft_gradphi_x)(i,j,k) = (*grid->ft_phi)(i,j,k)*idqx*n;
	  (*grid->ft_gradphi_y)(i,j,k) = (*grid->ft_phi)(i,j,k)*idqy*m;
	  (*grid->ft_gradphi_z)(i,j,k) = (*grid->ft_phi)(i,j,k)*idqz*l;
	  
	  
	}
      }
    }
    
    fftw_execute(grid->ps_grid->backward_gradphi_x);
    fftw_execute(grid->ps_grid->backward_gradphi_y);
    fftw_execute(grid->ps_grid->backward_gradphi_z);
    
    double factor = grid->boxgrid[0]*grid->boxgrid[1]*grid->boxgrid[2];
    (*grid->gradphi_x) /= factor;
    (*grid->gradphi_y) /= factor;
    (*grid->gradphi_z) /= factor;
  }


}


