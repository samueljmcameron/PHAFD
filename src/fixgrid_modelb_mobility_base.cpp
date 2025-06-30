#include <algorithm>

#include "utility.hpp"


#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "integrate.hpp"
#include "fixgrid_modelb_mobility_base.hpp"
#include "fixgrid_gradient.hpp"
#include "fixgrid_divergence.hpp"
#include "conjugate_noise_no_q.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;


FixGridModelBMobilityBase::FixGridModelBMobilityBase(PHAFD *phafd) : Fix(phafd), plan_set(false) {};



void FixGridModelBMobilityBase::init(const std::vector<std::string> &v_line)
{

  
  Fix::init(v_line);

  temp = -1;
  
  std::array<int,3> seeds;
  std::array<std::string,3> vlabels = {"noise_x","noise_y","noise_z"};

  int iarg = 1;
  for (; iarg < 4; iarg++)
    seeds.at(iarg-1)
      = utility::make_unique_seed(std::stoi(v_line.at(iarg)),
				  world,commbrick->me,commbrick->nprocs);



  while (iarg < v_line.size()) {

    if (v_line.at(iarg) == "temp") {
      temp = std::stod(v_line.at(iarg+1));
      iarg += 2;
    } else {
      throw std::runtime_error("Error: invalid fix grid/modelb/mobility command");
    }
  }

  if (temp < 0) 
    throw std::runtime_error("Error: invalid temp in fix grid/modelb/mobility command");


  std::vector<std::string> new_v_line;

  for (int i = 0; i < 3; i++) {
    conjugate_noise.at(i) = std::make_unique<ConjugateNoiseNoQ>(phafd);
    new_v_line.clear();
    new_v_line.push_back(vlabels.at(i));
    new_v_line.push_back("seed");
    new_v_line.push_back(std::to_string(seeds.at(i)));
    new_v_line.push_back("mobility");
    // set mobility to one since variable mobility takes care of
    // the actual prefactor.
    new_v_line.push_back("1.0"); 
    new_v_line.push_back("temp");
    new_v_line.push_back(std::to_string(temp));
    conjugate_noise.at(i)->readCoeffs(new_v_line);
  }

  

  for (int dim = 0; dim < 3; dim++)
    ft_rnoises[dim] = conjugate_noise.at(dim)->ft_array.get();
  
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

    for (int dim = 0; dim < 3; dim++)
      fftw_destroy_plan(backward_rnoises[dim]);

    fftw_destroy_plan(forward_mobility_deriv);
    fftw_destroy_plan(forward_sqrt_mobility);
    
  }

}

void FixGridModelBMobilityBase::setup()
{

  normalization = 1.0/(grid->ft_boxgrid[0]*grid->ft_boxgrid[1]*grid->ft_boxgrid[2]);
  int Nx = grid->boxgrid[0];
  int Ny = grid->boxgrid[1];
  int Nz = grid->boxgrid[2];

  inv_vol_element
    = Nx/domain->period[0]*Ny/domain->period[1]*Nz/domain->period[2]; 

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


    if (!rnoises[dim])
      rnoises[dim] = std::make_unique<fftwArr::array3D<double>>
	(world,"rnoises[" + std::to_string(dim) + "]" + name, Nx,Ny,Nz);
    
    backward_rnoises[dim]
      = fftw_mpi_plan_dft_c2r_3d(Nz,Ny,Nx,reinterpret_cast<fftw_complex*>
				 (ft_rnoises[dim]->data()),
				 rnoises[dim]->data(),
				 world,FFTW_MPI_TRANSPOSED_IN);

    
  }
  


  if (!sqrt_mobility)
    sqrt_mobility = std::make_unique<fftwArr::array3D<double>>
      (world,"sqrt_mobility" +  name,Nx,Ny,Nz);

  if (!ft_sqrt_mobility)
    ft_sqrt_mobility
      = std::make_unique<fftwArr::array3D<std::complex<double>>>
      (world,"ft_sqrt_mobility" +  name,Nx,Ny,Nz);

  
  forward_sqrt_mobility
    = fftw_mpi_plan_dft_r2c_3d(Nz,Ny,Nx,sqrt_mobility->data(),
				 reinterpret_cast<fftw_complex*>
				 (ft_sqrt_mobility->data()),
				 world, FFTW_MPI_TRANSPOSED_OUT);


  if (!mobility_deriv)
    mobility_deriv = std::make_unique<fftwArr::array3D<double>>
      (world,"mobility_deriv" +  name,Nx,Ny,Nz);


  if (!ft_mobility_deriv)
    ft_mobility_deriv
      = std::make_unique<fftwArr::array3D<std::complex<double>>>
      (world,"ft_mobility_deriv" +  name,Nx,Ny,Nz);

  
  forward_mobility_deriv
    = fftw_mpi_plan_dft_r2c_3d(Nz,Ny,Nx,mobility_deriv->data(),
				 reinterpret_cast<fftw_complex*>
				 (ft_mobility_deriv->data()),
				 world, FFTW_MPI_TRANSPOSED_OUT);

  plan_set = true;    
  
  
}

void FixGridModelBMobilityBase::reset_dt()
{

  dt = integrate->dt;
  for (int dim = 0; dim < 3; dim++)
    conjugate_noise.at(dim)->reset_dt(dt);
  
}


void FixGridModelBMobilityBase::start_of_step()
{

  // can compute these two things right away as they only require phi
  calculate_mobility_deriv();
  calculate_sqrt_mobility();


  // can compute the forward fourier transforms
  fftw_execute(grid->forward_phi);
  fftw_execute(forward_mobility_deriv);
  fftw_execute(forward_sqrt_mobility);

  

}


void FixGridModelBMobilityBase::compute_stochastic_drift()
{


  
  gradfix->calculate_gradient(ft_mobility_deriv.get());

  for (int i = 0; i < grid->chempot->Nz(); i++) 
    for (int j = 0; j < grid->chempot->Ny(); j++)
      for (int k = 0; k < grid->chempot->Nx(); k++) {
	
	for (int dim = 0; dim < 3; dim ++) // store stochastic drift term
	  (*flux[dim])(i,j,k)
	    = temp*(*gradfix->gradient[dim])(i,j,k)*inv_vol_element;

      }

}

void FixGridModelBMobilityBase::compute_usual_drift()
{
  double prefac;

  fftw_execute(grid->forward_chempot);

  gradfix->calculate_gradient(grid->ft_chempot.get());


  for (int i = 0; i < grid->chempot->Nz(); i++) 
    for (int j = 0; j < grid->chempot->Ny(); j++)
      for (int k = 0; k < grid->chempot->Nx(); k++) {
	
	prefac = (*sqrt_mobility)(i,j,k)*(*sqrt_mobility)(i,j,k);

	for (int dim = 0; dim < 3; dim ++) // store mobility*gradient
	  (*flux[dim])(i,j,k) += prefac*(*gradfix->gradient[dim])(i,j,k);


      }
}

void FixGridModelBMobilityBase::add_noise_to_phi()
{


  for (int dim = 0; dim < 3; dim++)
    conjugate_noise.at(dim)->update();

  divfix->calculate_divergence(ft_rnoises[0],ft_rnoises[1],
			       ft_rnoises[2]);

  // must do this after all ft_rnoises calculations are done.
  for (int dim = 0; dim < 3; dim++)
    fftw_execute(backward_rnoises[dim]);


  gradfix->calculate_gradient(ft_sqrt_mobility.get());
  
  int localNx = grid->phi->Nx();
  int localNy = grid->phi->Ny();
  int localNz = grid->phi->Nz();

  for (int i = 0; i < localNz; i++)
      for (int j = 0; j < localNy; j++)
	for (int k = 0; k < localNx; k++) {
	  for (int dim = 0; dim < 3; dim++)
	    (*grid->phi)(i,j,k)
	      += (*gradfix->gradient[dim])(i,j,k)*(*rnoises[dim])(i,j,k)*normalization;

	  (*grid->phi)(i,j,k)
	    += (*sqrt_mobility)(i,j,k)*(*divfix->divergence)(i,j,k);
	}  
}

void FixGridModelBMobilityBase::pre_final_integrate()
{

  // must calculate in this order, since stochastic drift resets flux
  compute_stochastic_drift();
  // whereas usual drift adds to flux
  compute_usual_drift();

  // fourier transform flux to then get divergence of it
  for (int dim = 0; dim < 3; dim++)
    fftw_execute(forward_flux[dim]);


  // calculate divergence of flux but leave it in fourier space
  divfix->calculate_divergence(ft_flux,false);
  
  return;
}


void FixGridModelBMobilityBase::final_integrate()
{

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

  add_noise_to_phi();
  
  return;
}





void FixGridModelBMobilityBase::point_update(int i , int j, int k)
{

  (*grid->ft_phi)(i,j,k)
    = ((*grid->ft_phi)(i,j,k)+dt*(*divfix->ft_divergence)(i,j,k))*normalization;
       // + (*grid->ft_noise)(i,j,k))*normalization;
  
  return;
  
}


void FixGridModelBMobilityBase::origin_update()
{

  (*grid->ft_phi)(0,0,0) *= normalization;

  return;
}

