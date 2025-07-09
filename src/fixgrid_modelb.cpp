#include <algorithm>

#include "utility.hpp"


#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "integrate.hpp"
#include "fixgrid_modelb.hpp"
#include "fixgrid_gradient.hpp"
#include "conjugate_noise.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;


FixGridModelB::FixGridModelB(PHAFD *phafd) : Fix(phafd) {};



void FixGridModelB::init(const std::vector<std::string> &v_line,
			 bool add_to_names)
/*
  v_line should have form
  fixname,seed

  followed by the key value pairs

  mobility, value
  temp, value

  in some order.
 */
{

  per_grid = true;
  Fix::init(v_line,add_to_names);

  double temp;
  mobility = temp = -1;


  // initialise conjugate class
  conjugate = std::make_unique<ConjugateNoise>(phafd);
  

  int seed = utility::make_unique_seed(std::stoi(v_line.at(1)),
				       world,commbrick->me,
				       commbrick->nprocs);


  std::vector<std::string> new_v_line;
  
  new_v_line.push_back("concentration");
  new_v_line.push_back("seed");
  new_v_line.push_back(std::to_string(seed));



  int iarg = 2;



  while (iarg < v_line.size()) {

    if (v_line.at(iarg) == "mobility") {
      new_v_line.push_back(v_line.at(iarg));
      new_v_line.push_back(v_line.at(iarg+1));
      mobility = std::stod(v_line.at(iarg+1));
      iarg += 2;
    } else if (v_line.at(iarg) == "temp") {
      new_v_line.push_back(v_line.at(iarg));
      new_v_line.push_back(v_line.at(iarg+1));
      temp = std::stod(v_line.at(iarg+1));
      iarg += 2;
    } else {
      throw std::runtime_error("Error: invalid fix grid/modelb command");
    }
  }


  if (mobility < 0) 
    throw std::runtime_error("Error: invalid mobility in fix grid/modelb command");

  if (temp < 0) 
    throw std::runtime_error("Error: invalid temp in fix grid/modelb command");

  
  conjugate->readCoeffs(new_v_line);
  
  complexFFTWarray.push_back(conjugate->ft_array.get());

  local0start = grid->phi->get_local0start();
  localNx = grid->phi->Nx();
  localNy = grid->phi->Ny();
  localNz = grid->phi->Nz();

  numberofcomponents=3;
  array.resize(localNx*localNy*localNz*numberofcomponents);

  new_v_line.clear();
  new_v_line.push_back(name+"_gradient");

  gradfix = std::make_unique<FixGridGradient>(phafd);
  gradfix->init(new_v_line,false);

  
  
}


void FixGridModelB::setup()
{

  gradfix->setup();
  normalization = 1.0/(grid->ft_boxgrid[0]*grid->ft_boxgrid[1]*grid->ft_boxgrid[2]);
  
}

void FixGridModelB::reset_dt()
{

  dt = integrate->dt;
  conjugate->reset_dt(dt);
  
}


void FixGridModelB::start_of_step()
{
  fftw_execute(grid->forward_phi);
}



void FixGridModelB::pre_final_integrate()
{


  fftw_execute(grid->forward_chempot);
  if (this_step) {
    gradfix->calculate_gradient(grid->ft_chempot.get());
    int count = 0;
    for (int i = 0; i < localNz; i++)
      for (int j = 0; j < localNy; j++)
	for (int k = 0; k < localNx; k++) {
	  array[count++] = -mobility*(*gradfix->gradient[0])(i,j,k);
	  array[count++] = -mobility*(*gradfix->gradient[1])(i,j,k);
	  array[count++] = -mobility*(*gradfix->gradient[2])(i,j,k);
	}
  }
    

  return;
}


void FixGridModelB::final_integrate()
{
  
  conjugate->update();
  int ft_local0start = grid->ft_phi->get_local0start();

  int ft_localNx = grid->ft_phi->Nx();
  int ft_localNy = grid->ft_phi->Ny();
  int ft_localNz = grid->ft_phi->Nz();



  
  if (commbrick->me == 0) {
    origin_update();

    for (int nx = 1; nx < ft_localNx; nx++)
      point_update(0,0,nx);
    
    for (int ny = 1; ny < ft_localNy; ny++)
      for (int nx = 0; nx < ft_localNx; nx++)
	point_update(0,ny,nx);

    
    for (int nz = 1; nz < ft_localNz; nz++)
      for (int ny = 0; ny < ft_localNy; ny++)
	for (int nx = 0; nx < ft_localNx; nx++)
	  point_update(nz,ny,nx);

  } else {
    for (int nz = 0; nz < ft_localNz; nz++)
      for (int ny = 0; ny < ft_localNy; ny++)
	for (int nx = 0; nx < ft_localNx; nx++)
	  point_update(nz,ny,nx);


  }
  
  return;
}


void FixGridModelB::post_final_integrate(bool invert_fft)
{
  if (invert_fft)

    fftw_execute(grid->backward_phi);
  
  return;
}





void FixGridModelB::point_update(int i , int j, int k)
{

  double qx,qy,qz,q2;
  
  qz = grid->qzs[i];
  qy = grid->qys[j];
  qx = domain->dqx()*k;

  q2 = qx*qx + qy*qy + qz*qz;

  (*grid->ft_phi)(i,j,k)
    = ((*grid->ft_phi)(i,j,k)-mobility*q2*dt*(*grid->ft_chempot)(i,j,k)
       + (*conjugate->ft_array)(i,j,k))*normalization;

  return;
  
}


void FixGridModelB::origin_update()
{

  (*grid->ft_phi)(0,0,0) *= normalization;

  return;
}

