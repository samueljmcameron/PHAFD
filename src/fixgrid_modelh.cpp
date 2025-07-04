#include <algorithm>

#include "utility.hpp"


#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "integrate.hpp"
#include "fixgrid_modelh.hpp"
#include "fixgrid_gradient.hpp"
#include "fixgrid_velocity.hpp"
#include "fixgrid_vdet.hpp"
#include "fixgrid_vtherm.hpp"
#include "conjugate_noise.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;


FixGridModelH::FixGridModelH(PHAFD *phafd) : Fix(phafd) {};



void FixGridModelH::init(const std::vector<std::string> &v_line)
/*
  v_line should have form
  
  fixname,conj_seed, vseedx, vseedy, vseedz
  
  followed by the key value pairs
  
  mobility, value
  temp, value
  viscosity, value
  
  in some order.
 */
{

  Fix::init(v_line);

  



  int conj_seed = utility::make_unique_seed(std::stoi(v_line.at(1)),
					    world,commbrick->me,
					    commbrick->nprocs);


  std::vector<int> velocity_seeds;
  int iarg = 2;
  for (; iarg < 5; iarg ++) 
    velocity_seeds.push_back(
			     utility::make_unique_seed(
				       std::stoi(v_line.at(iarg)
						 ),
				       world,commbrick->me,
				       commbrick->nprocs)
			     );

  mobility = -1;
  double temp = -1;
  double viscosity = -1;


  while (iarg < v_line.size()) {

    if (v_line[iarg] == "mobility") {
      mobility = std::stod(v_line[iarg+1]);
      iarg += 2;

    } else if (v_line[iarg] == "temp") {
      temp = std::stod(v_line[iarg+1]);
      iarg += 2;
    } else if (v_line[iarg] == "viscosity") {
      viscosity = std::stod(v_line[iarg+1]);
      iarg += 2;
    } else {
      throw std::runtime_error("Error: invalid fix grid/modelh command");
    }
  }


  if (mobility < 0) 
    throw std::runtime_error("Error: invalid mobility in fix grid/modelh command");

  if (temp < 0) 
    throw std::runtime_error("Error: invalid temp in fix grid/modelh command");

  if (viscosity < 0) 
    throw std::runtime_error("Error: invalid viscosity in fix grid/modelh command");
  

  
  conjugate = std::make_unique<ConjugateNoise>(phafd);

  std::vector<std::string> new_v_line;
  new_v_line.push_back("concentration");
  new_v_line.push_back("seed");
  new_v_line.push_back(std::to_string(conj_seed));

  new_v_line.push_back("mobility");
  new_v_line.push_back(std::to_string(mobility));
  new_v_line.push_back("temp");
  new_v_line.push_back(std::to_string(temp));  

  
  conjugate->readCoeffs(new_v_line);




  local_fixes.push_back(std::make_unique<FixGridVelocity>(phafd));
  
  // clear the array of strings to rewrite it for velocities
  new_v_line.clear();

  new_v_line.push_back(name + "_velocity");
  for (auto &seed : velocity_seeds) 
    new_v_line.push_back(std::to_string(seed));

  new_v_line.push_back(std::to_string(viscosity));
  new_v_line.push_back(std::to_string(temp));
  new_v_line.push_back("immediate");

  local_fixes.back()->init(new_v_line);

  local_fixes.push_back(std::make_unique<FixGridGradient>(phafd));

  // clear the array of strings to rewrite it for gradphi
  new_v_line.clear();
  new_v_line.push_back(name + "_gradphi");
  new_v_line.push_back("immediate");

  local_fixes.back()->init(new_v_line);

  /*

    fixgrid_v_dot_gradphi = std::make_unique<FixGridVdotGradPhi>(phafd);

  // clear the array of strings to rewrite it for v_dot_gradphi
  new_v_line.clear();

  new_v_line.push_back(name + "_v_dot_gradphi");
  for (auto &seed : velocity_seeds) 
    new_v_line.push_back(std::to_string(seed));

  new_v_line.push_back(std::to_string(viscosity));
  new_v_line.push_back(std::to_string(temp));
  new_v_line.push_back("immediate");

  
  fixgrid_v_dot_gradphi->init(new_v_line);  
  */


}


void FixGridModelH::setup()
{


  normalization = 1.0/(grid->ft_boxgrid[0]*grid->ft_boxgrid[1]*grid->ft_boxgrid[2]);

  for (auto &lf : local_fixes)
    lf->setup();

}

void FixGridModelH::reset_dt()
{

  dt = integrate->dt;
  conjugate->reset_dt(dt);
  for (auto &lf : local_fixes)
    lf->reset_dt();
}


void FixGridModelH::start_of_step()
{
  for (auto &lf : local_fixes)
    lf->start_of_step();


  fftw_execute(grid->forward_phi);
}

void FixGridModelH::initial_integrate()
{
  for (auto &lf : local_fixes)
    lf->initial_integrate();



}

void FixGridModelH::post_force()
{

  for (auto &lf : local_fixes)
    lf->post_force();

}


void FixGridModelH::pre_final_integrate()
{


  fftw_execute(grid->forward_chempot);
  for (auto &lf : local_fixes)
    lf->pre_final_integrate();


  //didnotintegrate = true;
  return;
}


void FixGridModelH::final_integrate()
{

  for (auto &lf : local_fixes)
    lf->final_integrate();
  
  conjugate->update();

  local0start = grid->ft_phi->get_local0start();

  localNx = grid->ft_phi->Nx();
  localNy = grid->ft_phi->Ny();
  localNz = grid->ft_phi->Nz();

  if (commbrick->me == 0) {
    origin_update();

    for (int nx = 1; nx < localNx; nx++)
      point_update(0,0,nx);
    
    for (int ny = 1; ny < localNy; ny++)
      for (int nx = 0; nx < localNx; nx++)
	point_update(0,ny,nx);

    
    for (int nz = 1; nz < localNz; nz++)
      for (int ny = 0; ny < localNy; ny++)
	for (int nx = 0; nx < localNx; nx++)
	  point_update(nz,ny,nx);

  } else {
    for (int nz = 0; nz < localNz; nz++)
      for (int ny = 0; ny < localNy; ny++)
	for (int nx = 0; nx < localNx; nx++)
	  point_update(nz,ny,nx);


  }
  
  //didnotintegrate = false;
  
  return;
}


void FixGridModelH::post_final_integrate(bool invert_fft)
{


  for (auto &lf : local_fixes)
    lf->post_final_integrate(invert_fft);

  if (invert_fft) {
    fftw_execute(grid->backward_phi);


    // CAREFUL HERE!! NEED TO ENSURE LOCAL are correct if outputting to array
    localNx = grid->phi->Nx();
    localNy = grid->phi->Ny();
    localNz = grid->phi->Nz();



    for (int i = 0; i < localNz; i++)
      for (int j = 0; j < localNy; j++)
	for (int k = 0; k < localNx; k++) 
	  for (int dim = 0; dim < 3; dim++)
	    (*grid->phi)(i,j,k) -=
	      (*grid->velocity[dim])(i,j,k)*(*grid->gradphi[dim])(i,j,k)*dt;
    
  }

  
  return;
}


void FixGridModelH::end_of_step()
{
  for (auto &lf : local_fixes)
    lf->end_of_step();



}


void FixGridModelH::point_update(int i , int j, int k)
{

  double qx,qy,qz,q2;
  
  qz = grid->qzs[i];
  qy = grid->qys[j];
  qx = domain->dqx()*k;

  q2 = qx*qx + qy*qy + qz*qz;


  // partially update the fourier space phi, ignoring those pieces
  // which are better off being in real space update.

  // for example, not include v.grad(phi) part here because it is
  // already present in real space.
  (*grid->ft_phi)(i,j,k)
    = ((*grid->ft_phi)(i,j,k)-mobility*q2*dt*(*grid->ft_chempot)(i,j,k)
       + (*conjugate->ft_array)(i,j,k))*normalization;
  
  return;
  
}


void FixGridModelH::origin_update()
{

  (*grid->ft_phi)(0,0,0) *= normalization;

  return;
}

