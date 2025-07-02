#include <algorithm>

#include "utility.hpp"


#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "integrate.hpp"
#include "fixgrid_vdet.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;


FixGridVdet::FixGridVdet(PHAFD *phafd) : Fix(phafd) {};



void FixGridVdet::init(const std::vector<std::string> &v_line)
/*
  v_line should take form:
  fixname,viscosity
 */
{

  Fix::init(v_line);

  viscosity = std::stod(v_line.at(1));

  if (v_line.size()==3) {
    if (v_line.at(2) == "immediate") immediate_ifft = true;
    else
      throw std::runtime_error("invalid argument in fix/vdet");
  }

  vdet_x = grid->vdet[0].get();
  vdet_y = grid->vdet[1].get();
  vdet_z = grid->vdet[2].get();
  ft_vdet_x = grid->ft_vdet[0].get();
  ft_vdet_y = grid->ft_vdet[1].get();
  ft_vdet_z = grid->ft_vdet[2].get();

  
  
}

  
void FixGridVdet::setup()
{

  grid->set_qs(*ft_vdet_x);
  normalization = (grid->ft_boxgrid[0]*grid->ft_boxgrid[1]*grid->ft_boxgrid[2]);

  
}



void FixGridVdet::pre_final_integrate()
{

  // store chemical potential x gradphi
  for (int i = 0; i < grid->chempot->Nz(); i++) 
    for (int j = 0; j < grid->chempot->Ny(); j++) 
      for (int k = 0; k < grid->chempot->Nx(); k++) {
	(*vdet_x)(i,j,k) = (*grid->gradphi[0])(i,j,k)*(*grid->chempot)(i,j,k);
	(*vdet_y)(i,j,k) = (*grid->gradphi[1])(i,j,k)*(*grid->chempot)(i,j,k);
	(*vdet_z)(i,j,k) = (*grid->gradphi[2])(i,j,k)*(*grid->chempot)(i,j,k);
      }

  // then do a fft to compute the fourier space version

  fftw_execute(grid->forward_vdet[0]);
  fftw_execute(grid->forward_vdet[1]);
  fftw_execute(grid->forward_vdet[2]);

  // then multiply fourier space version by the Oseen tensor

  int localNx = ft_vdet_x->Nx();
  int localNy = ft_vdet_x->Ny();
  int localNz = ft_vdet_x->Nz();

  for (int nz = 0; nz < localNz; nz ++)
    for (int ny = 0; ny < localNy; ny ++)
      for (int nx = 0; nx < localNx; nx++)
	set_vdet(nz,ny,nx);


  if (immediate_ifft) {
    for (int i = 0; i < 3; i++) 
      fftw_execute(grid->backward_vdet[i]);
  }


}


void FixGridVdet::post_final_integrate(bool invert_fft) {
  if (immediate_ifft) return;


  /* IMPORTANT!!!! THE TWO OPERATIONS DONE BELOW MUST BE DONE
     IN SEPARATE LOOPS, SINCE THE fftw_execute COMMAND SWAPS
     vdet[1] AND vdet[2] !!!
  */

  if (invert_fft) {
    for (int i = 0; i < 3; i++) 
      fftw_execute(grid->backward_vdet[i]);
    
    /* DO NOT COMBINE THIS LOOP WITH THE ABOVE LOOP!!! */
    for (int i = 0; i < 3; i++) {
      (*grid->vdet[i]) /= normalization;
    }
  }

  
}


// setting the values of vdet in fourier space, transposes are accounted for here.

void FixGridVdet::set_vdet(int i, int j, int k) {

  double qx,qy,qz,q2,Txx,Txy,Txz,Tyy,Tyz,Tzz;

  std::complex<double> tmp_x,tmp_y,tmp_z;


  qy = grid->qzs[i];
  qz = grid->qys[j];
  qx = domain->dqx()*k;


  q2 = qx*qx + qy*qy + qz*qz;

  double denom = viscosity;
  if (immediate_ifft)
    denom *= normalization;


  if (q2 == 0) {
    (*ft_vdet_x)(i,j,k) = 0.0;
    (*ft_vdet_y)(i,j,k) = 0.0;
    (*ft_vdet_z)(i,j,k) = 0.0;


  } else {



    Txx = (1.0-qx*qx/q2)/(q2*denom);
    Txy = (-qx*qy/q2)/(q2*denom);
    Txz = (-qx*qz/q2)/(q2*denom);
    Tyy = (1.0-qy*qy/q2)/(q2*denom);
    Tyz = (-qy*qz/q2)/(q2*denom);
    Tzz = (1.0-qz*qz/q2)/(q2*denom);

    tmp_x = (*ft_vdet_x)(i,j,k);
    tmp_y = (*ft_vdet_y)(i,j,k);
    tmp_z = (*ft_vdet_z)(i,j,k);

    (*ft_vdet_x)(i,j,k) = (Txx*tmp_x+Txy*tmp_y
			   + Txz*tmp_z);
    (*ft_vdet_y)(i,j,k) = (Txy*tmp_x+Tyy*tmp_y
			   + Tyz*tmp_z);
    (*ft_vdet_z)(i,j,k) = (Txz*tmp_x+Tyz*tmp_y
			   + Tzz*tmp_z);
    

  }

  
  return ;

}
