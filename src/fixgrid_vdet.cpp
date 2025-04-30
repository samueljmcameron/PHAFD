#include <algorithm>

#include "utility.hpp"


#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "integrate.hpp"
#include "fixgrid_vdet.hpp"
#include "conjugate_noise.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;


FixGridVdet::FixGridVdet(PHAFD *phafd) : Fix(phafd) {};



void FixGridVdet::init(const std::vector<std::string> &v_line)
/*
  v_line should take form:
  fixname,seedx,seedy,seedz,viscosity,temperature
 */
{

  Fix::init(v_line);

  std::vector<std::string> new_v_line;

  // build string vector for concentration

  std::array<std::string,3> vlabels = {"vx","vy","vz"};

  std::array<int,3> seeds;

  int iarg = 1;
  for (; iarg < 4; iarg++)
    seeds.at(iarg-1) = utility::make_unique_seed(std::stoi(v_line.at(iarg)),
						 world,commbrick->me,commbrick->nprocs);

  viscosity = std::stod(v_line.at(4));
  temp = std::stod(v_line.at(5));
  
  for (int i = 0; i < 3; i++) {
    conjugate_vnoise.at(i) = std::make_unique<ConjugateNoise>(phafd);
    new_v_line.clear();
    new_v_line.push_back(vlabels.at(i));
    new_v_line.push_back("seed");
    new_v_line.push_back(std::to_string(seeds.at(i)));
    new_v_line.push_back("viscosity");
    new_v_line.push_back(std::to_string(viscosity));
    new_v_line.push_back("temp");
    new_v_line.push_back(std::to_string(temp));

    
    conjugate_vnoise.at(i)->readCoeffs(new_v_line);
  }


  viscosity = conjugate_vnoise.at(0)->damping;
  
  ft_Znoise_x = grid->ft_Znoise[0].get();
  ft_Znoise_y = grid->ft_Znoise[1].get();
  ft_Znoise_z = grid->ft_Znoise[2].get();

  ft_vdet_x = grid->ft_vdet[0].get();
  ft_vdet_y = grid->ft_vdet[1].get();
  ft_vdet_z = grid->ft_vdet[2].get();
  
}

  
void FixGridVdet::setup()
{
  conjugate_vnoise.at(0)->copy_qs(qys,qzs);

}

void FixGridVdet::reset_dt()
{

  dt = integrate->dt;

  for (auto &cv : conjugate_vnoise)
    cv->reset_dt(dt);
  
}


void FixGridVdet::start_of_step()
{
  for (auto &cv : conjugate_vnoise)
    cv->update();


  // and compute v_therm in fourier space;
  compute_vdet();
  // and inverse fourier transform to get vdet in real space
  for (int i = 0; i < 3; i++) 
    fftw_execute(grid->backward_vdet[i]);

  
}



void FixGridVdet::compute_vdet() {

  int localNx = ft_Znoise_x->Nx();
  int localNy = ft_Znoise_x->Ny();
  int localNz = ft_Znoise_x->Nz();
  
  for (int nz = 0; nz < localNz; nz++) 
    for (int ny = 0; ny < localNy; ny++) 
      for (int nx = 0; nx < localNx; nx++)
	set_vdet(nz,ny,nx);

}


// setting the values of vdet in fourier space, transposes are accounted for here.

void FixGridVdet::set_vdet(int i, int j, int k) {

  double qx,qy,qz,q2,Txx,Txy,Txz,Tyy,Tyz,Tzz;


  qz = qzs[i];
  qy = qys[j];
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

    (*ft_vdet_x)(i,j,k) = (Txx*(*ft_Znoise_x)(i,j,k)+Txz*(*ft_Znoise_y)(i,j,k)
			      + Txy*(*ft_Znoise_z)(i,j,k))/dt;
    (*ft_vdet_y)(i,j,k) = (Txz*(*ft_Znoise_x)(i,j,k)+Tzz*(*ft_Znoise_y)(i,j,k)
			      + Tyz*(*ft_Znoise_z)(i,j,k))/dt;
    (*ft_vdet_z)(i,j,k) = (Txy*(*ft_Znoise_x)(i,j,k)+Tyz*(*ft_Znoise_y)(i,j,k)
			      + Tyy*(*ft_Znoise_z)(i,j,k))/dt;


  }
  
  return ;

}
