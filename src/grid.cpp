#include <iostream>
#include <random>
#include "grid.hpp"
#include "comm_brick.hpp"

#include "utility.hpp"
#include "domain.hpp"
#include "read_dump.hpp"
#include "fftw_arr/array3d.hpp"


using namespace PHAFD_NS;

Grid::Grid(PHAFD *phafd) : Pointers(phafd) {

  /* set for any grid */
  gridset = false;
  gridpopulated = false;

  


  /* concentration arrays */

  phi_set = false;
  
  phi = nullptr;
  chempot = nullptr;
  ft_phi = nullptr;
  ft_chempot = nullptr;
  for (int i = 0; i < 3; i++) {
    gradphi[i] = nullptr;
    ft_gradphi[i] = nullptr;
  }


  /* velocity arrays */

  velocity_set = false;
  
  for (int i = 0; i < 3; i++) {
    velocity[i] = nullptr;
    ft_velocity[i] = nullptr;
    vtherm[i] = nullptr;
    ft_vtherm[i] = nullptr;
  }


  /* model H arrays */

  modelH_set = false;

  vtherm_dot_gradphi = nullptr;
  ft_vtherm_dot_gradphi = nullptr;
  for (int i = 0; i < 3; i++) {
    ft_Znoise[i] = nullptr;
  }
  
}


Grid::~Grid() {

  if (phi_set) {
    fftw_destroy_plan(forward_phi);
    fftw_destroy_plan(backward_phi);
    fftw_destroy_plan(forward_chempot);
    fftw_destroy_plan(backward_chempot);
    for (int i = 0; i < 3; i++) {
      fftw_destroy_plan(backward_gradphi[i]);
    }

  }

  if (velocity_set) {
    for (int i = 0; i < 3; i++) {
      fftw_destroy_plan(forward_velocity[i]);
      fftw_destroy_plan(backward_velocity[i]);
      fftw_destroy_plan(forward_vtherm[i]);
      fftw_destroy_plan(backward_vtherm[i]);

      
    }
  }

}

void Grid::create(const std::vector<std::string> &v_line)
{

  int Nx = std::stoi(v_line.at(0));
  int Ny = std::stoi(v_line.at(1));
  int Nz = std::stoi(v_line.at(2));

  boxgrid = {Nx,Ny,Nz};
  // since using transposed fftw methods, need ft_boxgrid to swap Nz and Ny
  ft_boxgrid = {Nx,Nz,Ny};

  int iarg = 3;
  
  while (iarg < v_line.size()) {
    
    if (v_line.at(iarg) == "concentration") {
      create_concentration(Nx,Ny,Nz);
      iarg += 1;
    } else {
      throw std::invalid_argument("Incompatible style type in grid_style.");
    }
  }
  
  gridset = true;
  
}

void Grid::populate(const std::vector<std::string> &v_line)
{

  std::vector<std::string> new_v_line = v_line;

  if (new_v_line.at(0) == "constant") {

    if (new_v_line.at(1) == "concentration") {
      
      double average = std::stod(new_v_line.at(2));
      double variance = std::stod(new_v_line.at(3));

      int seed = std::stoi(new_v_line.at(4));
      seed = utility::make_unique_seed(seed,world,commbrick->me,commbrick->nprocs);

      noisy_constant(phi.get(),average,variance,seed);
      
    } else {
      throw std::runtime_error("Invalid options for grid_populate.");
    }
    
  } else if (new_v_line.at(0) == "read_dump") {
      
    new_v_line.at(0) = "grid";

    ReadDump read_dump(phafd);

    read_dump.init(new_v_line);
    read_dump.process_attributes();
    
  } else {
    throw std::runtime_error("Invalid options for grid_populate.");
  }
  


  gridpopulated=true;
}


double Grid::dx() const
{
  return domain->period[0]/boxgrid[0];
}

double Grid::dy() const
{
  return domain->period[1]/boxgrid[1];
}


double Grid::dz() const
{
  return domain->period[2]/boxgrid[2];
}



void Grid::create_concentration(int Nx, int Ny, int Nz)
{
  phi_set = true;
  
  std::array<std::string,3> listxyz = {"x","y","z"};

  if (!phi) 
    phi = std::make_unique<fftwArr::array3D<double>
			   >(world,"concentration",Nx,Ny,Nz);
  if (!ft_phi)
    ft_phi = std::make_unique<fftwArr::array3D<std::complex<double>>
			      >(world,"ft_concentration",Nx,Nz,Ny);

  if (!chempot) 
    chempot = std::make_unique<fftwArr::array3D<double>
				 >(world,"chemPot",Nx,Ny,Nz);
  if (!ft_chempot)
    ft_chempot = std::make_unique<fftwArr::array3D<std::complex<double>>
				    >(world,"ft_chemPot",Nx,Nz,Ny);

  
  for (int i = 0; i < 3; i++) {
    if (!gradphi[i])
      gradphi[i] = std::make_unique<fftwArr::array3D<double>
				    >(world,"gradphi_"+listxyz[i],Nx,Ny,Nz);

    if (!ft_gradphi[i])
      ft_gradphi[i] = std::make_unique<fftwArr::array3D<std::complex<double>>
				       >(world,"ft_gradphi_"+listxyz[i],Nx,Nz,Ny);

  }


  
  forward_phi = fftw_mpi_plan_dft_r2c_3d(Nz,Ny,Nx,phi->data(),
					 reinterpret_cast<fftw_complex*>
					 (ft_phi->data()),
					 world, FFTW_MPI_TRANSPOSED_OUT);
  
  backward_phi = fftw_mpi_plan_dft_c2r_3d(Nz,Ny,Nx,reinterpret_cast<fftw_complex*>
					  (ft_phi->data()), phi->data(),
					  world,FFTW_MPI_TRANSPOSED_IN);
  
  forward_chempot = fftw_mpi_plan_dft_r2c_3d(Nz,Ny,Nx,
					       chempot->data(),
					       reinterpret_cast<fftw_complex*>
					       (ft_chempot->data()),
					       world, FFTW_MPI_TRANSPOSED_OUT);
  
  backward_chempot = fftw_mpi_plan_dft_c2r_3d(Nz,Ny,Nx,
						reinterpret_cast<fftw_complex*>
						(ft_chempot->data()),
						chempot->data(),world,
						FFTW_MPI_TRANSPOSED_IN);



  backward_gradphi[0] = fftw_mpi_plan_dft_c2r_3d(Nz,Ny,Nx,
						 reinterpret_cast<fftw_complex*>
						 (ft_gradphi[0]->data()),
						 gradphi[0]->data(),world,
						 FFTW_MPI_TRANSPOSED_IN);


  // NOTE HERE THE SWAP IN Z AND Y! THIS IS NOT A BUG!!!
  backward_gradphi[1] = fftw_mpi_plan_dft_c2r_3d(Nz,Ny,Nx,
						 reinterpret_cast<fftw_complex*>
						 (ft_gradphi[1]->data()),
						 gradphi[2]->data(),world, // <---HERE!! NOT A BUG!
						 FFTW_MPI_TRANSPOSED_IN);
  // NOTE HERE THE SWAP IN Z AND Y! THIS IS NOT A BUG!!!
  backward_gradphi[2] = fftw_mpi_plan_dft_c2r_3d(Nz,Ny,Nx,
						 reinterpret_cast<fftw_complex*>
						 (ft_gradphi[2]->data()),
						 gradphi[1]->data(),world,// <---HERE!! NOT A BUG!
						 FFTW_MPI_TRANSPOSED_IN);


  
  return;
  
}


void Grid::noisy_constant(fftwArr::array3D<double> * array,
			  double average, double variance, int seed)
{

  if (!array)
    throw std::runtime_error("Cannot create constant array.");
  std::uniform_real_distribution<double> real_dist(-0.5,0.5);
    
  std::mt19937 gen;
    
  gen.seed(seed);
	  
  for (int i = 0; i < array->Nz(); i++) 
    for (int j = 0; j < array->Ny(); j++) 
      for (int k = 0; k < array->Nx(); k++) 
	(*array)(i,j,k) = average + variance*real_dist(gen);
  
  return;
}
