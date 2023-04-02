
#include "conjugate_volfrac.hpp"

#include <cmath>
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;

ConjugateVolFrac::ConjugateVolFrac(PHAFD *phafd)
  : Conjugate(phafd),real_dist(-0.5,0.5)
{

  ft_phi = grid->ft_phi.get();
  ft_chempot = grid->ft_chempot.get();
  
  if (!ft_phi || !ft_chempot)
    throw std::runtime_error("Calling conjugate/volfrac on grid that "
			     "doesn't have concentration grid_style.");

  
  seed_flag = false;
  ft_array = ft_phi;
  
  setup();
  
}


void ConjugateVolFrac::readCoeffs(const std::vector<std::string> &v_line)
{


  mobility = temp = volFH = gamma = 1;

  int iarg = 0;

  while (iarg < v_line.size()) {

    if (v_line[iarg] == "mobility") {
      mobility = std::stod(v_line[iarg+1]);
      iarg += 2;

    } else if (v_line[iarg] == "temp") {
      temp = std::stod(v_line[iarg+1]);
      iarg += 2;
    } else if (v_line[iarg] == "volFH") {
      volFH = std::stod(v_line[iarg+1]);
      iarg += 2;
    } else if (v_line[iarg] == "gamma") {
      gamma = std::stod(v_line[iarg+1]);
      iarg += 2;
    } else if (v_line[iarg] == "seed") {
      seed = std::stod(v_line[iarg+1]);
      iarg += 2;
      seed_flag = true;
      
    } else {
      throw std::runtime_error("Error: invalid conjugate/volfrac command");
    }
  }

  if (!seed_flag)
    throw std::runtime_error("Error: need seed specified in conjugate/volfrac command");

  gen.seed(seed);

  
  

}

void ConjugateVolFrac::reset_dt(double timestep)
{
  dt = timestep;
  double invLcubed = 1.0/(domain->period[0]*domain->period[1]*domain->period[2]);

  normalization = 1.0/(grid->ft_boxgrid[0]*grid->ft_boxgrid[1]*grid->ft_boxgrid[2]);


  complexprefactor = sqrt(12*temp*mobility*invLcubed);
  realprefactor = sqrt(24*temp*mobility*invLcubed);

  sqrtdt = sqrt(dt);


}


void ConjugateVolFrac::complex_update(int i , int j, int k)
{

  double qx,qy,qz,q2;
  
  qz = qzs[i];
  qy = qys[j];
  qx = domain->dqx()*k;

  q2 = qx*qx + qy*qy + qz*qz;




  noise.real(complexprefactor*sqrt(q2)*sqrtdt*real_dist(gen));
  noise.imag(complexprefactor*sqrt(q2)*sqrtdt*real_dist(gen));

  (*ft_phi)(i,j,k)
    = ((*ft_phi)(i,j,k)-mobility*q2*dt*((*ft_chempot)(i,j,k)
					+temp/volFH*gamma*q2*(*ft_phi)(i,j,k))
       )*normalization + noise;
  
  return;
  
}


void ConjugateVolFrac::real_update(int i, int j, int k)
{


  double qx,qy,qz,q2;

  qz = qzs[i];
  qy = qys[j];
  qx = domain->dqx()*k;
  
  q2 = qx*qx + qy*qy + qz*qz;

  noise = realprefactor*sqrt(q2)*sqrtdt*real_dist(gen);
  
  (*ft_phi)(i,j,k)
    = ((*ft_phi)(i,j,k)-mobility*q2*dt*((*ft_chempot)(i,j,k)
					+temp/volFH*gamma*q2*(*ft_phi)(i,j,k))
       )*normalization + noise;
  
  return;
}

void ConjugateVolFrac::origin_update()
{

  (*ft_phi)(0,0,0) *= normalization;

  return;
}

