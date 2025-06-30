
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

  damping = temp = -1;



  name = v_line.at(0);
  
  int iarg = 1;

  
  while (iarg < v_line.size()) {

    if (v_line.at(iarg) == "mobility" || v_line.at(iarg) == "viscosity") {
      damping = std::stod(v_line.at(iarg+1));
      iarg += 2;

    } else if (v_line.at(iarg) == "temp") {
      temp = std::stod(v_line.at(iarg+1));
      iarg += 2;
    } else if (v_line.at(iarg) == "seed") {
      seed = std::stod(v_line.at(iarg+1));
      iarg += 2;
      seed_flag = true;
      
    } else {
      throw std::runtime_error("Error: invalid conjugate/noise command");
    }
  }

  if (damping < 0)
    throw std::runtime_error("Error: invalid damping in conjugate/noise command");

  if (temp < 0)
    throw std::runtime_error("Error: invalid damping in conjugate/noise command");
  

  gen.seed(seed);

  setup();  
  

}

void ConjugateNoise::reset_dt(double timestep)
{
  dt = timestep;
  double invLcubed = 1.0/(domain->period[0]*domain->period[1]
			  *domain->period[2]);

  double un_normalization = (grid->ft_boxgrid[0]*grid->ft_boxgrid[1]
			     *grid->ft_boxgrid[2]);


  complexprefactor = sqrt(12*temp*damping*invLcubed)*un_normalization;
  realprefactor = sqrt(24*temp*damping*invLcubed)*un_normalization;
  
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

