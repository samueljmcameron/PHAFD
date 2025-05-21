#include <algorithm>

#include "utility.hpp"

#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "integrate.hpp"
#include "fixgrid_gradphi.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;


FixGridGradPhi::FixGridGradPhi(PHAFD *phafd) : Fix(phafd) {
  ftphi_flag = false;
};



void FixGridGradPhi::init(const std::vector<std::string> &v_line)
{

  Fix::init(v_line);

  for (int iarg = 1; iarg < v_line.size(); iarg++) {
    if (v_line.at(iarg) == "ftphi")  ftphi_flag = true;
    else if (v_line.at(iarg) == "immediate") immediate_ifft = true;
    else
      throw std::runtime_error("invalid fix/gradphi");
  }
  
}

  
void FixGridGradPhi::setup()
{
  normalization = (grid->boxgrid[0]*grid->boxgrid[1]*grid->boxgrid[2]);

}


void FixGridGradPhi::start_of_step()
{

  if (ftphi_flag == true) {
    fftw_execute(grid->forward_phi);
  }

}

void FixGridGradPhi::post_force() {
  
  const int local0start = grid->ft_phi->get_local0start();
  const int globalNy = grid->ft_boxgrid[1];
  const int globalNz = grid->ft_boxgrid[2];
  
  
  std::complex<double> idqx(0,domain->dqx());
  std::complex<double> idqy(0,domain->dqy());
  std::complex<double> idqz(0,domain->dqz());
  
  
  if (immediate_ifft) {
    idqx /= normalization;
    idqy /= normalization;
    idqz /= normalization;
  }
  
  
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
	
	(*grid->ft_gradphi[0])(i,j,k) = (*grid->ft_phi)(i,j,k)*idqx*n;
	(*grid->ft_gradphi[1])(i,j,k) = (*grid->ft_phi)(i,j,k)*idqy*m;
	(*grid->ft_gradphi[2])(i,j,k) = (*grid->ft_phi)(i,j,k)*idqz*l;
	
      }
    }
  }
  
  if (immediate_ifft) {
    for (int i = 0; i < 3; i++) 
      fftw_execute(grid->backward_gradphi[i]);
  }
  
  
}



void FixGridGradPhi::post_final_integrate() {

  if (immediate_ifft) return;

  /* IMPORTANT!!!! THE TWO OPERATIONS DONE BELOW MUST BE DONE
     IN SEPARATE LOOPS, SINCE THE fftw_execute COMMAND SWAPS
     gradphi[1] AND gradphi[2] !!!
  */

  for (int i = 0; i < 3; i++)
    fftw_execute(grid->backward_gradphi[i]);


  /* DO NOT COMBINE THIS LOOP WITH THE ABOVE LOOP!!! */
  for (int i = 0; i < 3; i++) {
    (*grid->gradphi[i]) /= normalization;
  }
  
}
