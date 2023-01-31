#include <iostream>

#include "grid.hpp"
#include "comm_brick.hpp"

#include "utility.hpp"
#include "domain.hpp"
#include "ps_pde/grid.hpp"

using namespace PHAFD_NS;

Grid::Grid(PHAFD *phafd) : Pointers(phafd) {

  ps_grid = nullptr;
  boxgrid = nullptr;
  ft_boxgrid = nullptr;
  phi = nullptr;
  nonlinear = nullptr;
  ft_phi = nullptr;
  ft_nonlinear = nullptr;
  ft_noise = nullptr;

  
}

Grid::~Grid() = default;


void Grid::create(const std::vector<std::string> &v_line)
{
  ps_grid = std::make_unique<psPDE::Grid>(v_line,world);

  boxgrid = ps_grid->boxgrid.data();
  ft_boxgrid = ps_grid->ft_boxgrid.data();

  phi = ps_grid->phi.get();
  nonlinear = ps_grid->nonlinear.get();

  ft_phi = ps_grid->ft_phi.get();
  ft_nonlinear = ps_grid->ft_nonlinear.get();
  ft_noise = ps_grid->ft_noise.get();

  
}

/* need to call after commbrick->setup() */
void Grid::populate(const std::vector<std::string> &v_line)
{

  std::vector<std::string> new_v_line = v_line;


  if (new_v_line[0] == "constant" && new_v_line[1] == "concentration") {
    int seed = std::stoi(new_v_line.back());


    seed = utility::make_unique_seed(seed,world,commbrick->me,commbrick->nprocs);

    new_v_line.back() = std::to_string(seed);
    
  } else if (new_v_line[0] != "read")
    
    throw std::runtime_error("Invalid options for grid_populate.");
  
  
  ps_grid->populate(new_v_line);
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
