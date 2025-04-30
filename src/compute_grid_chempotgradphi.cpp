#include <algorithm>

#include "utility.hpp"


#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "integrate.hpp"
#include "compute_grid_chempotgradphi.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;


ComputeGridVdet::ComputeGridVdet(PHAFD *phafd) : Compute(phafd) {};



void ComputeGridVdet::init(const std::vector<std::string> &v_line)
/*
  v_line should take form:
  fixname,seedx,seedy,seedz,viscosity,temperature
 */
{

  Compute::init(v_line);

  std::vector<std::string> new_v_line;

  // build string vector for concentration

  viscosity = std::stod(v_line.at(1));
  
}

  
void ComputeGridVdet::setup()
{
  //conjugate_vnoise.at(0)->copy_qs(qys,qzs);

}

void ComputeGridVdet::reset_dt()
{
  
}


void ComputeGridVdet::start_of_step()
{
  
}



void ComputeGridVdet::compute_vdet() {

  int localNx = ft_Znoise_x->Nx();
  int localNy = ft_Znoise_x->Ny();
  int localNz = ft_Znoise_x->Nz();
  
  for (int nz = 0; nz < localNz; nz++) 
    for (int ny = 0; ny < localNy; ny++) 
      for (int nx = 0; nx < localNx; nx++)
	set_vdet(nz,ny,nx);

}


// setting the values of vdet in fourier space, transposes are accounted for here.

void ComputeGridVdet::set_vdet(int i, int j, int k) {

  double qx,qy,qz,q2,Txx,Txy,Txz,Tyy,Tyz,Tzz;


  qz = grid->qzs[i];
  qy = grid->qys[j];
  qx = domain->dqx()*k;
  
  q2 = qx*qx + qy*qy + qz*qz;


  if (q2 == 0) {
    (*ft_vdet_x)(i,j,k) = 0.0;
    (*ft_vdet_y)(i,j,k) = 0.0;
    (*ft_vdet_z)(i,j,k) = 0.0;
  } else {



    Txx = (1.0-qx*qx/q2)/(q2*viscosity);
    Txy = (-qx*qy/q2)/(q2*viscosity);  
    Txz = (-qx*qz/q2)/(q2*viscosity);
    Tyy = (1-qy*qy/q2)/(q2*viscosity);
    Tyz = (-qy*qz/q2)/(q2*viscosity);
    Tzz = (1-qz*qz/q2)/(q2*viscosity);


    // need to do a swap here (qy <-> qz) since computing transposed
    // fourier functions, the below LOOKS LIKE IT HAS BUGS BUT IT DOES NOT!!

    (*ft_chempot_gradphi_x)(i,j,k) = (Txx*(*ft_Znoise_x)(i,j,k)+Txz*(*ft_Znoise_y)(i,j,k)
			      + Txy*(*ft_Znoise_z)(i,j,k))/dt;
    (*ft_chempot_gradphi_y)(i,j,k) = (Txz*(*ft_Znoise_x)(i,j,k)+Tzz*(*ft_Znoise_y)(i,j,k)
			      + Tyz*(*ft_Znoise_z)(i,j,k))/dt;
    (*ft_chempot_gradphi_z)(i,j,k) = (Txy*(*ft_Znoise_x)(i,j,k)+Tyz*(*ft_Znoise_y)(i,j,k)
			      + Tyy*(*ft_Znoise_z)(i,j,k))/dt;


  }
  
  return ;

}
