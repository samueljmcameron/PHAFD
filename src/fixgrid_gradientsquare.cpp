#include <algorithm>
#include <iostream>
#include "utility.hpp"


#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "fixgrid_laplacian.hpp"
#include "fixgrid_gradientsquare.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;


FixGridGradientSquare::FixGridGradientSquare(PHAFD *phafd) : Fix(phafd) {

};



void FixGridGradientSquare::init(const std::vector<std::string> &v_line,
				 bool add_to_names)
{

  Fix::init(v_line,add_to_names);
  prefac = std::stod(v_line.at(1));

  ft_phi_flag = false;
  if (v_line.size() == 3) {
    if (v_line.at(2) == "ft_phi")
      ft_phi_flag = true;
    else
      throw  std::runtime_error("Fix grid/gradientsquare error.");
  } else if (v_line.size() > 3)
    throw std::runtime_error("Fix grid/gradientsquare error.");

  std::vector<std::string> new_v_line;
  new_v_line.push_back(name+"_laplacian");

  laplacian = std::make_unique<FixGridLaplacian>(phafd);
  laplacian->init(new_v_line,false);

  
}

void FixGridGradientSquare::setup() {

  laplacian->setup();
  
}

void FixGridGradientSquare::post_force()
{

  if (ft_phi_flag) {
    fftw_execute(grid->forward_phi);
  }
  laplacian->calculate_laplacian(grid->ft_phi.get());

  
  for (int i = 0; i < grid->chempot->Nz(); i++) 
    for (int j = 0; j < grid->chempot->Ny(); j++)
      for (int k = 0; k < grid->chempot->Nx(); k++) 
	(*grid->chempot)(i,j,k) 
	  -= prefac*(*laplacian->laplacian)(i,j,k);

}

