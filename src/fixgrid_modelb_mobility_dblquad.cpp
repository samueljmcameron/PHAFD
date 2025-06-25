#include <algorithm>

#include "utility.hpp"

#include "fixgrid_modelb_mobility_dblquad.hpp"


#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "integrate.hpp"
#include "conjugate_noise.hpp"
#include "fixgrid_gradient.hpp"
#include "fixgrid_divergence.hpp"
#include "fftw_arr/array3d.hpp"


using namespace PHAFD_NS;


FixGridModelBMobilityDBLquad::FixGridModelBMobilityDBLquad(PHAFD *phafd) : FixGridModelBMobilityBase(phafd) {};



void FixGridModelBMobilityDBLquad::init(const std::vector<std::string> &v_line)
/*
  v_line should take form:
  fixname,prefactor
 */
{

  if (v_line.end()[-2] != "mobility_prefactor")
    throw std::runtime_error("Error: invalid fix grid/modelb/mobility/dblquad command");
  
  prefactor = std::stod(v_line.back());
  std::vector<std::string> new_v_line = v_line;
  new_v_line.pop_back();
  new_v_line.pop_back();
  FixGridModelBMobilityBase::init(new_v_line);

  
}


void FixGridModelBMobilityDBLquad::calculate_sqrt_mobility()
{
  for (int i = 0; i < sqrt_mobility->Nz(); i++) 
    for (int j = 0; j < sqrt_mobility->Ny(); j++)
      for (int k = 0; k < sqrt_mobility->Nx(); k++)
	(*sqrt_mobility)(i,j,k)
	  = sqrt(prefactor)*std::abs((*grid->phi)(i,j,k)*(1-(*grid->phi)(i,j,k)));

}
