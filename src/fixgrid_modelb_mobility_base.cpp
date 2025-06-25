#include <algorithm>

#include "utility.hpp"


#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "integrate.hpp"
#include "fixgrid_modelb_mobility_base.hpp"
#include "fixgrid_gradient.hpp"
#include "fixgrid_divergence.hpp"
#include "conjugate_noise.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;


FixGridModelBMobilityBase::FixGridModelBMobilityBase(PHAFD *phafd) : Fix(phafd), plan_set(false) {};



void FixGridModelBMobilityBase::init(const std::vector<std::string> &v_line)
/*
  v_line should have form
  fixname,seed

  followed by the key value pairs

  mobility, value
  temp, value

  in some order.
 */
{

  
  Fix::init(v_line);

  double temp = -1;


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

    if (v_line.at(iarg) == "temp") {
      new_v_line.push_back(v_line.at(iarg));
      new_v_line.push_back(v_line.at(iarg+1));
      temp = std::stod(v_line.at(iarg+1));
      iarg += 2;
    } else {
      throw std::runtime_error("Error: invalid fix grid/modelb/mobility command");
    }
  }

  if (temp < 0) 
    throw std::runtime_error("Error: invalid temp in fix grid/modelb/mobility command");

  new_v_line.push_back("mobility");
  new_v_line.push_back("1.0");

  
  conjugate->readCoeffs(new_v_line);



  new_v_line.clear();
  new_v_line.push_back(name+"_gradient");

  gradfix = std::make_unique<FixGridGradient>(phafd);
  gradfix->init(new_v_line);



  new_v_line.clear();
  new_v_line.push_back(name+"_divergence");

  divfix = std::make_unique<FixGridDivergence>(phafd);
  divfix->init(new_v_line);


  
}


FixGridModelBMobilityBase::~FixGridModelBMobilityBase()
{
  if (plan_set) {
    for (int dim = 0; dim < 3; dim++)
      fftw_destroy_plan(forward_flux[dim]);
  }

}

void FixGridModelBMobilityBase::setup()
{

  normalization = 1.0/(grid->ft_boxgrid[0]*grid->ft_boxgrid[1]*grid->ft_boxgrid[2]);
  int Nx = grid->boxgrid[0];
  int Ny = grid->boxgrid[1];
  int Nz = grid->boxgrid[2];

  gradfix->setup();
  divfix->setup();

  for (int dim = 0; dim < 3; dim++) {
    if (!flux[dim])
      flux[dim] = std::make_unique<fftwArr::array3D<double>>
	(world,"flux[" + std::to_string(dim) + "]" +  name,Nx,Ny,Nz);
    if (!ft_flux[dim])
      ft_flux[dim]
	= std::make_unique<fftwArr::array3D<std::complex<double>>>
	(world,"ft_flux[" + std::to_string(dim) + "]" +  name,Nx,Ny,Nz);
    
    forward_flux[dim]
      = fftw_mpi_plan_dft_r2c_3d(Nz,Ny,Nx,flux[dim]->data(),
				 reinterpret_cast<fftw_complex*>
				 (ft_flux[dim]->data()),
				 world, FFTW_MPI_TRANSPOSED_OUT);
    
    if (!grad_sqrt_mobility[dim])
      grad_sqrt_mobility[dim]
	= std::make_unique<fftwArr::array3D<double>>
	(world,"grad_sqrt_mobility["
	 + std::to_string(dim) + "]" +  name,Nx,Ny,Nz);
    
    
    
  }
  
  plan_set = true;  
  

  if (!sqrt_mobility)
    sqrt_mobility = std::make_unique<fftwArr::array3D<double>>
      (world,"sqrt_mobility" +  name,Nx,Ny,Nz);

  if (!mobility_deriv)
    mobility_deriv = std::make_unique<fftwArr::array3D<double>>
      (world,"mobility_deriv" +  name,Nx,Ny,Nz);  
  
  
}

void FixGridModelBMobilityBase::reset_dt()
{

  dt = integrate->dt;
  conjugate->reset_dt(dt);
  
}


void FixGridModelBMobilityBase::start_of_step()
{
  fftw_execute(grid->forward_phi);
}



void FixGridModelBMobilityBase::pre_final_integrate()
{


  fftw_execute(grid->forward_chempot);

  gradfix->calculate_gradient(grid->ft_chempot.get());

  double prefac;

  calculate_sqrt_mobility();

  for (int i = 0; i < grid->chempot->Nz(); i++) 
    for (int j = 0; j < grid->chempot->Ny(); j++)
      for (int k = 0; k < grid->chempot->Nx(); k++) {
	
	prefac = (*sqrt_mobility)(i,j,k)*(*sqrt_mobility)(i,j,k);

	for (int dim = 0; dim < 3; dim ++) // store mobility*gradient
	  (*flux[dim])(i,j,k) = prefac*(*gradfix->gradient[dim])(i,j,k);


      }
  // fourier transform flux to then get divergence of it
  for (int dim = 0; dim < 3; dim++)
    fftw_execute(forward_flux[dim]);



  // calculate divergence of flux but leave it in fourier space
  divfix->calculate_divergence(ft_flux,false);
  
  return;
}


void FixGridModelBMobilityBase::final_integrate()
{
  
  conjugate->update();

  int local0start = grid->ft_phi->get_local0start();

  int localNx = grid->ft_phi->Nx();
  int localNy = grid->ft_phi->Ny();
  int localNz = grid->ft_phi->Nz();
  
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
  
  return;
}


void FixGridModelBMobilityBase::post_final_integrate()
{

  fftw_execute(grid->backward_phi);
  
  return;
}





void FixGridModelBMobilityBase::point_update(int i , int j, int k)
{

  (*grid->ft_phi)(i,j,k)
    = ((*grid->ft_phi)(i,j,k)+dt*(*divfix->ft_divergence)(i,j,k)
       + (*grid->ft_noise)(i,j,k))*normalization;
  
  return;
  
}


void FixGridModelBMobilityBase::origin_update()
{

  (*grid->ft_phi)(0,0,0) *= normalization;

  return;
}

