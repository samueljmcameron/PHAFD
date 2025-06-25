#include <algorithm>
#include <iostream>
#include "utility.hpp"

#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "integrate.hpp"
#include "fixgrid_laplacianphi.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;


FixGridLaplacianPhi::FixGridLaplacianPhi(PHAFD *phafd) : Fix(phafd) {
  ftphi_flag = false;
  
};



void FixGridLaplacianPhi::init(const std::vector<std::string> &v_line)
{

  Fix::init(v_line);

  for (int iarg = 1; iarg < v_line.size(); iarg++) {
    if (v_line.at(iarg) == "ftphi")  ftphi_flag = true;
    else if (v_line.at(iarg) == "immediate") immediate_ifft = true;
    else
      throw std::runtime_error("invalid fix/laplacianphi");
  }

  
}

  
void FixGridLaplacianPhi::setup()
{
  normalization = (grid->boxgrid[0]*grid->boxgrid[1]*grid->boxgrid[2]);

}

void FixGridLaplacianPhi::start_of_step() {
  if (ftphi_flag == true) {
    fftw_execute(grid->forward_phi);
  }

}

void FixGridLaplacianPhi::post_final_integrate() {
  if (immediate_ifft) return;
  fftw_execute(grid->backward_laplacianphi);
  (*grid->laplacianphi) /= normalization;

}

void FixGridLaplacianPhi::post_force()
{

  const int local0start = grid->ft_phi->get_local0start();
  const int globalNy = grid->ft_boxgrid[1];
  const int globalNz = grid->ft_boxgrid[2];
  
  double dqy = domain->dqy();
  double dqz = domain->dqz();
  double dqx = domain->dqx();
  
  if (immediate_ifft) {
    dqx /= sqrt(normalization);
    dqy /= sqrt(normalization);
    dqz /= sqrt(normalization);
  }
  
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
	
	neg_q2 = -(dqx*n*dqx*n + dqz*m*dqz*m + dqy*l*dqy*l);
	
	
	(*grid->ft_laplacianphi)(i,j,k) = (*grid->ft_phi)(i,j,k)*neg_q2;
	
      }
    }
  }
  
  if (immediate_ifft)
    fftw_execute(grid->backward_laplacianphi);
  
}


