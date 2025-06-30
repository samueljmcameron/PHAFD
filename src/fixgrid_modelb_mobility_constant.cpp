#include <algorithm>

#include "utility.hpp"

#include "fixgrid_modelb_mobility_constant.hpp"


#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "integrate.hpp"
#include "conjugate_noise.hpp"
#include "fixgrid_gradient.hpp"
#include "fixgrid_divergence.hpp"
#include "fftw_arr/array3d.hpp"


using namespace PHAFD_NS;


FixGridModelBMobilityConstant::FixGridModelBMobilityConstant(PHAFD *phafd) : FixGridModelBMobilityBase(phafd) {};



void FixGridModelBMobilityConstant::init(const std::vector<std::string> &v_line)
/*
  v_line should have form
  fixname,seedx,seedy,seedz

  followed by the key value pair

  temp, value

  followed by the key value pair

  mobility_prefactor, value

 */
{

  if (v_line.end()[-2] != "mobility_prefactor")
    throw std::runtime_error("Error: invalid fix grid/modelb/mobility/constant command");
  
  prefactor = std::stod(v_line.back());
  std::vector<std::string> new_v_line = v_line;
  new_v_line.pop_back();
  new_v_line.pop_back();
  FixGridModelBMobilityBase::init(new_v_line);

  
}


void FixGridModelBMobilityConstant::calculate_sqrt_mobility()
{
  for (int i = 0; i < sqrt_mobility->Nz(); i++) 
    for (int j = 0; j < sqrt_mobility->Ny(); j++)
      for (int k = 0; k < sqrt_mobility->Nx(); k++)
	(*sqrt_mobility)(i,j,k)
	  = sqrt(prefactor);

}


void FixGridModelBMobilityConstant::calculate_mobility_deriv()
{
  for (int i = 0; i < mobility_deriv->Nz(); i++) 
    for (int j = 0; j < mobility_deriv->Ny(); j++)
      for (int k = 0; k < mobility_deriv->Nx(); k++)
	(*mobility_deriv)(i,j,k)
	  = 0.0;

}
