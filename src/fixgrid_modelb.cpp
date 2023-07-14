#include <algorithm>

#include "utility.hpp"


#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "integrate.hpp"
#include "fixgrid_modelb.hpp"
#include "fixgrid_gradphi.hpp"
#include "conjugate_volfrac.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;


FixGridModelB::FixGridModelB(PHAFD *phafd) : Fix(phafd) {};



void FixGridModelB::init(const std::vector<std::string> &v_line)
/*
  v_line should have form
  fixname,seed,mobility,temp,volFH,gamma
 */
{

  Fix::init(v_line);

  conjugate = std::make_unique<ConjugateVolFrac>(phafd);
  
  std::vector<std::string> new_v_line;

  new_v_line.push_back(v_line.at(1));

  for (int i = 2; i < 6; i++)
    new_v_line.push_back(v_line.at(i));

  conjugate->readCoeffs(new_v_line);


  // re-write new_vline to include a fix for computing the gradient of phi
  new_v_line.clear();

  new_v_line.push_back(v_line.at(0)+std::string("gradphi"));
  new_v_line.push_back("every");
  fixgridgradphi = std::make_unique<FixGridGradPhi>(phafd);
  fixgridgradphi->init(new_v_line);


}


void FixGridModelB::setup()
{

}

void FixGridModelB::reset_dt()
{

  dt = integrate->dt;
  conjugate->reset_dt(dt);
  
}


void FixGridModelB::start_of_step()
{
  fixgridgradphi->start_of_step();
}



void FixGridModelB::pre_final_integrate()
{


  fftw_execute(grid->forward_chempot);

  didnotintegrate = true;
  return;
}


void FixGridModelB::final_integrate()
{


  conjugate->update();
  didnotintegrate = false;
  
  return;
}


void FixGridModelB::post_final_integrate()
{

  fftw_execute(grid->backward_phi);

  if (didnotintegrate) {
    double factor = grid->boxgrid[0]*grid->boxgrid[1]*grid->boxgrid[2];
    
    (*grid->phi) /= factor;
  }
  
  return;
}



