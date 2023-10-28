#include <algorithm>

#include "utility.hpp"


#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "integrate.hpp"
#include "fixgrid_modelb.hpp"
#include "fixgrid_gradphi.hpp"
#include "conjugate_noise.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;


FixGridModelB::FixGridModelB(PHAFD *phafd) : Fix(phafd) {};



void FixGridModelB::init(const std::vector<std::string> &v_line)
/*
  v_line should have form
  fixname,seed,mobility,temp,volFH,gamma
 */
{

  ft_phi = grid->ft_phi.get();
  ft_chempot = grid->ft_chempot.get();
  
  Fix::init(v_line);

  conjugate = std::make_unique<ConjugateNoise>(phafd);
  
  std::vector<std::string> new_v_line;


  int seed = std::stoi(v_line.at(1));
  seed = utility::make_unique_seed(seed,world,commbrick->me,commbrick->nprocs);

  new_v_line.push_back("concentration");
  new_v_line.push_back("seed");
  new_v_line.push_back(std::to_string(seed));


  mobility = temp = volFH = gamma = 1;




  std::vector<std::string> another_new_line;

  int iarg = 2;

  
  while (iarg < v_line.size()) {

    if (v_line[iarg] == "mobility") {
      mobility = std::stod(v_line[iarg+1]);
      new_v_line.push_back(v_line[iarg]);
      new_v_line.push_back(v_line[iarg+1]);
      iarg += 2;

    } else if (v_line[iarg] == "temp") {
      temp = std::stod(v_line[iarg+1]);
      new_v_line.push_back(v_line[iarg]);
      new_v_line.push_back(v_line[iarg+1]);
      
      iarg += 2;
    } else if (v_line[iarg] == "volFH") {
      volFH = std::stod(v_line[iarg+1]);
      iarg += 2;
    } else if (v_line[iarg] == "gamma") {
      gamma = std::stod(v_line[iarg+1]);
      iarg += 2;
    } else {
      throw std::runtime_error("Error: invalid fix grid/modelb command");
    }
  }

  
  
  conjugate->readCoeffs(new_v_line);




  

  // re-write new_vline to include a fix for computing the gradient of phi
  new_v_line.clear();

  new_v_line.push_back(v_line.at(0)+std::string("gradphi"));
  new_v_line.push_back("every");
  fixgridgradphi = std::make_unique<FixGridGradPhi>(phafd);
  fixgridgradphi->init(new_v_line);


}


void FixGridModelB::setup()
{

  conjugate->copy_qs(qys,qzs);
  normalization = 1.0/(grid->ft_boxgrid[0]*grid->ft_boxgrid[1]*grid->ft_boxgrid[2]);
  
}

void FixGridModelB::reset_dt()
{

  dt = integrate->dt;
  conjugate->reset_dt(dt);
  
}


void FixGridModelB::start_of_step()
{
  fixgridgradphi->start_of_step();
}



void FixGridModelB::pre_final_integrate()
{


  fftw_execute(grid->forward_chempot);

  didnotintegrate = true;
  return;
}


void FixGridModelB::final_integrate()
{


  
  conjugate->update();

  int local0start = ft_phi->get_local0start();

  if (commbrick->me == 0) {
    origin_update();

    for (int nx = 1; nx < ft_phi->Nx(); nx++)
      point_update(0,0,nx);
    
    for (int ny = 1; ny < ft_phi->Ny(); ny++)
      for (int nx = 0; nx < ft_phi->Nx(); nx++)
	point_update(0,ny,nx);

    
    for (int nz = 1; nz < ft_phi->Nz(); nz++)
      for (int ny = 0; ny < ft_phi->Ny(); ny++)
	for (int nx = 0; nx < ft_phi->Nx(); nx++)
	  point_update(nz,ny,nx);

  } else {
    for (int nz = 0; nz < ft_phi->Nz(); nz++)
      for (int ny = 0; ny < ft_phi->Ny(); ny++)
	for (int nx = 0; nx < ft_phi->Nx(); nx++)
	  point_update(nz,ny,nx);


  }
  
  didnotintegrate = false;
  
  return;
}


void FixGridModelB::post_final_integrate()
{

  fftw_execute(grid->backward_phi);

  if (didnotintegrate) {
    double factor = grid->boxgrid[0]*grid->boxgrid[1]*grid->boxgrid[2];
    
    (*grid->phi) /= factor;
  }
  
  return;
}





void FixGridModelB::point_update(int i , int j, int k)
{

  double qx,qy,qz,q2;
  
  qz = qzs[i];
  qy = qys[j];
  qx = domain->dqx()*k;

  q2 = qx*qx + qy*qy + qz*qz;

  (*ft_phi)(i,j,k)
    = ((*ft_phi)(i,j,k)-mobility*q2*dt*((*ft_chempot)(i,j,k)
					+temp/volFH*gamma*q2*(*ft_phi)(i,j,k))
       )*normalization + (*grid->ft_noise)(i,j,k);
  
  return;
  
}


void FixGridModelB::origin_update()
{

  (*ft_phi)(0,0,0) *= normalization;

  return;
}

