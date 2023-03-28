#include <iostream>

#include "grid.hpp"
#include "comm_brick.hpp"

#include "utility.hpp"
#include "domain.hpp"
#include "ps_pde/grid.hpp"

using namespace PHAFD_NS;

Grid::Grid(PHAFD *phafd) : Pointers(phafd) {
  gridset = false;
  gridpopulated = false;
  ps_grid = nullptr;
  boxgrid = nullptr;
  ft_boxgrid = nullptr;
  phi = nullptr;
  chempot = nullptr;
  gradphi_x = nullptr;
  gradphi_y = nullptr;
  gradphi_z = nullptr;
  ft_phi = nullptr;
  ft_chempot = nullptr;
  ft_noise = nullptr;
  ft_gradphi_x = nullptr;
  ft_gradphi_y = nullptr;
  ft_gradphi_z = nullptr;
  
}

Grid::~Grid() = default;


void Grid::create(const std::vector<std::string> &v_line)
{
  ps_grid = std::make_unique<psPDE::Grid>(v_line,world);

  boxgrid = ps_grid->boxgrid.data();
  ft_boxgrid = ps_grid->ft_boxgrid.data();

  phi = ps_grid->phi.get();
  chempot = ps_grid->chempot.get();
  gradphi_x = ps_grid->gradphi_x.get();
  gradphi_y = ps_grid->gradphi_y.get();
  gradphi_z = ps_grid->gradphi_z.get();

  ft_phi = ps_grid->ft_phi.get();
  ft_chempot = ps_grid->ft_chempot.get();
  ft_noise = ps_grid->ft_noise.get();

  ft_gradphi_x = ps_grid->ft_gradphi_x.get();
  ft_gradphi_y = ps_grid->ft_gradphi_y.get();
  ft_gradphi_z = ps_grid->ft_gradphi_z.get();
  
  gridset = true;
  
}

void Grid::populate(const std::vector<std::string> &v_line)
{

  std::vector<std::string> new_v_line = v_line;

  if  (new_v_line.at(0) == "sinx" || new_v_line.at(0) == "siny"
       || new_v_line.at(0) == "sinz") {

    sinusoidal_grid(new_v_line.at(0),std::stod(new_v_line.at(1)));

    
  } else if (new_v_line.at(0) == "constant" && new_v_line.at(1) == "concentration") {
    int seed = std::stoi(new_v_line.back());


    seed = utility::make_unique_seed(seed,world,commbrick->me,commbrick->nprocs);

    new_v_line.back() = std::to_string(seed);

    ps_grid->populate(new_v_line);
    
  } else if (new_v_line.at(0) == "read") {
    ps_grid->populate(new_v_line);
  } else
    
    throw std::runtime_error("Invalid options for grid_populate.");
  
  


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


void Grid::sinusoidal_grid(const std::string & stype,double period)
{

  if (!phi)
    throw std::runtime_error("Cannot create sinusoidal phi.");


  if (stype == "sinx") {

    double Ox = domain->boxlo[0];

    
    for (int i = 0; i < phi->Nz(); i++)  {
      for (int j = 0; j < phi->Ny(); j++) {
	for (int k = 0; k < phi->Nx(); k++) { 
	  (*phi)(i,j,k) = sin(period*(k*dx()+Ox));
	}
      }
    }
    
  } else if (stype == "siny") {
    double y;
    double Oy = domain->boxlo[1];

    for (int i = 0; i < phi->Nz(); i++)  {
      for (int j = 0; j < phi->Ny(); j++) {
	y = j*dy()+Oy;
	for (int k = 0; k < phi->Nx(); k++)  {
	  (*phi)(i,j,k) = sin(period*y);
	}
      }
    }

  } else if (stype == "sinz") {

    int local0start = phi->get_local0start();
    double z;
    double Oz = domain->boxlo[2];

    for (int i = 0; i < phi->Nz(); i++) {
      z = (i+local0start)*dz()+Oz;
      for (int j = 0; j < phi->Ny(); j++)  {
	for (int k = 0; k < phi->Nx(); k++) {
	  (*phi)(i,j,k) = sin(period*z);
	}
      }
    }
  } 

  return;
				      
}
