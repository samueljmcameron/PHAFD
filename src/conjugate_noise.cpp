
#include "conjugate_noise.hpp"
#include "grid.hpp"
#include <cmath>
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;

ConjugateNoise::ConjugateNoise(PHAFD *phafd)
  : Conjugate(phafd),real_dist(-0.5,0.5)
{
  seed_flag = false;
  
}


void ConjugateNoise::readCoeffs(const std::vector<std::string> &v_line)
{

  damping = temp = 1;

  
  if (v_line.at(0) == "concentration") {
    ft_array = grid->ft_noise.get();
  } else if (v_line.at(0) == "vx") {
    ft_array = grid->ft_Znoise[0].get();
  } else if (v_line.at(0) == "vy") {
    ft_array = grid->ft_Znoise[1].get();
  } else if (v_line.at(0) == "vz") {
    ft_array = grid->ft_Znoise[2].get();
  } else
    throw std::runtime_error("Invalid keyword in conjugate/noise.");



  int iarg = 1;

  
  while (iarg < v_line.size()) {

    if (v_line[iarg] == "mobility" || v_line[iarg] == "viscosity") {
      damping = std::stod(v_line[iarg+1]);
      iarg += 2;

    } else if (v_line[iarg] == "temp") {
      temp = std::stod(v_line[iarg+1]);
      iarg += 2;
    } else if (v_line[iarg] == "seed") {
      seed = std::stod(v_line[iarg+1]);
      iarg += 2;
      seed_flag = true;
      
    } else {
      throw std::runtime_error("Error: invalid conjugate/noise command");
    }
  }


  gen.seed(seed);

  setup();  
  

}

void ConjugateNoise::reset_dt(double timestep)
{
  dt = timestep;
  double invLcubed = 1.0/(domain->period[0]*domain->period[1]*domain->period[2]);

  normalization = 1.0/(grid->ft_boxgrid[0]*grid->ft_boxgrid[1]*grid->ft_boxgrid[2]);


  complexprefactor = sqrt(12*temp*damping*invLcubed);
  realprefactor = sqrt(24*temp*damping*invLcubed);
  
  sqrtdt = sqrt(dt);


}


void ConjugateNoise::complex_update(int i , int j, int k)
{

  double qx,qy,qz,q2;
  
  qz = grid->qzs[i];
  qy = grid->qys[j];
  qx = domain->dqx()*k;

  q2 = qx*qx + qy*qy + qz*qz;



  (*ft_array)(i,j,k).real(complexprefactor*sqrt(q2)*sqrtdt*real_dist(gen));
  (*ft_array)(i,j,k).imag(complexprefactor*sqrt(q2)*sqrtdt*real_dist(gen));

  
  return;
  
}


void ConjugateNoise::real_update(int i, int j, int k)
{


  double qx,qy,qz,q2;

  qz = grid->qzs[i];
  qy = grid->qys[j];
  qx = domain->dqx()*k;
  
  q2 = qx*qx + qy*qy + qz*qz;

  (*ft_array)(i,j,k).real(realprefactor*sqrt(q2)*sqrtdt*real_dist(gen));
  
  return;
}

void ConjugateNoise::origin_update()
{

  (*ft_array)(0,0,0) = 0;

  return;
}

