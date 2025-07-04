#include <algorithm>

#include "utility.hpp"


#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "integrate.hpp"
#include "fixgrid_vtherm.hpp"
#include "conjugate_noise.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;


FixGridVtherm::FixGridVtherm(PHAFD *phafd) : Fix(phafd) {};



void FixGridVtherm::init(const std::vector<std::string> &v_line)
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

  if (v_line.size() == 7) {
    if (v_line.at(6) == "immediate")
      immediate_ifft = true;
    else
      throw std::runtime_error("invalid argument in fix/vtherm");
  }
  
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
  
  ft_Znoise_x = conjugate_vnoise.at(0)->ft_array.get();
  ft_Znoise_y = conjugate_vnoise.at(1)->ft_array.get();
  ft_Znoise_z = conjugate_vnoise.at(2)->ft_array.get();

  ft_vtherm_x = grid->ft_vtherm[0].get();
  ft_vtherm_y = grid->ft_vtherm[1].get();
  ft_vtherm_z = grid->ft_vtherm[2].get();
  
}

  
void FixGridVtherm::setup()
{
  normalization = (grid->boxgrid[0]*grid->boxgrid[1]*grid->boxgrid[2]);

}

void FixGridVtherm::reset_dt()
{

  dt = integrate->dt;

  for (auto &cv : conjugate_vnoise)
    cv->reset_dt(dt);
  
}


void FixGridVtherm::start_of_step()
{
  for (auto &cv : conjugate_vnoise)
    cv->update();


  // and compute v_therm in fourier space;
  compute_vtherm();



}


void FixGridVtherm::post_final_integrate(bool invert_fft) {

  if (immediate_ifft) return;
  if (invert_fft) 
    for (int i = 0; i < 3; i++) {
      fftw_execute(grid->backward_vtherm[i]);
      (*grid->vtherm[i]) /= normalization;
    }
}



void FixGridVtherm::compute_vtherm() {

  localNx = ft_Znoise_x->Nx();
  localNy = ft_Znoise_x->Ny();
  localNz = ft_Znoise_x->Nz();
  
  for (int nz = 0; nz < localNz; nz++) 
    for (int ny = 0; ny < localNy; ny++) 
      for (int nx = 0; nx < localNx; nx++)
	set_vtherm(nz,ny,nx);

  if (immediate_ifft)
    for (int i = 0; i < 3; i++) 
      fftw_execute(grid->backward_vtherm[i]);

  
}


// setting the values of vtherm in fourier space, transposes are accounted for here.

void FixGridVtherm::set_vtherm(int i, int j, int k) {

  double qx,qy,qz,q2,Txx,Txy,Txz,Tyy,Tyz,Tzz;


  qz = grid->qzs[i];
  qy = grid->qys[j];
  qx = domain->dqx()*k;
  
  q2 = qx*qx + qy*qy + qz*qz;


  double denom = viscosity;
  if (immediate_ifft)
    denom *= normalization;

  
  if (q2 == 0) {
    (*ft_vtherm_x)(i,j,k) = 0.0;
    (*ft_vtherm_y)(i,j,k) = 0.0;
    (*ft_vtherm_z)(i,j,k) = 0.0;
  } else {



    Txx = (1.0-qx*qx/q2)/(q2*denom);
    Txy = (-qx*qy/q2)/(q2*denom);  
    Txz = (-qx*qz/q2)/(q2*denom);
    Tyy = (1-qy*qy/q2)/(q2*denom);
    Tyz = (-qy*qz/q2)/(q2*denom);
    Tzz = (1-qz*qz/q2)/(q2*denom);


    (*ft_vtherm_x)(i,j,k) = (Txx*(*ft_Znoise_x)(i,j,k)+Txy*(*ft_Znoise_y)(i,j,k)
			      + Txz*(*ft_Znoise_z)(i,j,k))/dt;
    (*ft_vtherm_y)(i,j,k) = (Txy*(*ft_Znoise_x)(i,j,k)+Tyy*(*ft_Znoise_y)(i,j,k)
			      + Tyz*(*ft_Znoise_z)(i,j,k))/dt;
    (*ft_vtherm_z)(i,j,k) = (Txz*(*ft_Znoise_x)(i,j,k)+Tyz*(*ft_Znoise_y)(i,j,k)
			      + Tzz*(*ft_Znoise_z)(i,j,k))/dt;


  }
  
  return ;

}
