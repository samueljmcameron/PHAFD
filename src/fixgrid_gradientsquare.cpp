#include <algorithm>

#include "utility.hpp"


#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"

#include "fixgrid_gradientsquare.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;


FixGridGradientSquare::FixGridGradientSquare(PHAFD *phafd) : FixGridLaplacianPhi(phafd) {

};



void FixGridGradientSquare::init(const std::vector<std::string> &v_line)
{


  FixGridLaplacianPhi::init(v_line);

  prefac = std::stod(v_line.at(1));

}



void FixGridGradientSquare::post_force()
{

  FixGridLaplacianPhi::post_force();
  for (int i = 0; i < grid->chempot->Nz(); i++) 
    for (int j = 0; j < grid->chempot->Ny(); j++)
      for (int k = 0; k < grid->chempot->Nx(); k++) 
	(*grid->chempot)(i,j,k) 
	  -= prefac*(*grid->laplacianphi)(i,j,k);

}

