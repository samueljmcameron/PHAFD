#include <algorithm>

#include "utility.hpp"


#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "integrate.hpp"
#include "fixgrid_modelh.hpp"
#include "conjugate_volfrac.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;


FixGridModelH::FixGridModelH(PHAFD *phafd) : Fix(phafd) {};



void FixGridModelH::init(const std::vector<std::string> &v_line)
/*
  v_line should take form:

  fixname,seedc,seedx,seedy,seedz,mobility,viscosity,temperature

 */
{

  Fix::init(v_line);

  std::vector<std::string> new_v_line;


  int iarg = 1;
  for (; iarg < 5; iarg++)
    seeds.at(iarg-1) = std::stoi(v_line.at(iarg));


  seeds.at(0) = utility::make_unique_seed(seeds.at(0),world,commbrick->me,commbrick->nprocs);
  

  mobility = temp = volFH = gamma = viscosity = 1;

  while (iarg < v_line.size()) {

    if (v_line[iarg] == "mobility") {
      mobility = std::stod(v_line[iarg+1]);
      iarg += 2;
    } else if (v_line[iarg] == "temp") {
      temp = std::stod(v_line[iarg+1]);
      
      iarg += 2;
    } else if (v_line[iarg] == "volFH") {
      volFH = std::stod(v_line[iarg+1]);
      iarg += 2;
    } else if (v_line[iarg] == "gamma") {
      gamma = std::stod(v_line[iarg+1]);
      iarg += 2;
    } else if (v_line[iarg] == "viscosity") {
      gamma = std::stod(v_line[iarg+1]);
      iarg += 2;
    } else {
      throw std::runtime_error("Error: invalid fix grid/modelb command");
    }
  }


  conjugate_phinoise = std::make_unique<ConjugateNoise>(phafd);
  // build string vector for concentration noise
  new_v_line.push_back("concentration");
  new_v_line.push_back("seed");
  new_v_line.push_back(std::to_string(seeds.at(0)));
  new_v_line.push_back("mobility");
  new_v_line.push_back(std::to_string(mobility));
  new_v_line.push_back("temp");
  new_v_line.push_back(std::to_string(temp));

  conjugate_phinoise->readCoeffs(new_v_line);


  new_v_line.clear();
  new_v_line.push_back(v_line.at(0) + std::string("_vtherm"));
  for (int i = 1; i < 4; i++)
    new_v_line.push_back(std::to_string(seeds.at(i)));
  new_v_line.push_back("viscosity");
  new_v_line.push_back(std::to_string(viscosity));
  new_v_line.push_back("temp");
  new_v_line.push_back(std::to_string(temp));

  fixgridvtherm = std::make_unique<FixGridVtherm>(phafd);
  fixgridvtherm->init(new_v_line);
  
  // re-write new_vline to include a fix for computing the gradient of phi
  new_v_line.clear();

  new_v_line.push_back(v_line.at(0)+std::string("_gradphi"));
  new_v_line.push_back("every");
  fixgridgradphi = std::make_unique<FixGridGradPhi>(phafd);
  fixgridgradphi->init(new_v_line);

  
}

  
void FixGridModelH::setup()
{

  fixgridvtherm->setup();
  normalization = 1.0/(grid->ft_boxgrid[0]*grid->ft_boxgrid[1]*grid->ft_boxgrid[2]);


}

void FixGridModelH::reset_dt()
{

  dt = integrate->dt;


  conjugate_phinoise->reset_dt(dt);

  
  fixgridvtherm->reset_dt();
}


void FixGridModelH::start_of_step()
{
  // compute grad phi in real space, also phi in fourier space
  fixgridgradphi->start_of_step();

  // compute noises in fourier space
  conjugate_phinoise->update();

  // and compute v_therm in fourier space;
  fixgridvtherm->start_of_step();

  // then compute vtherm.gradphi in real space

  compute_vtherm_dot_gradphi();
  
  // and the grad\tilde{phi} in real space
  compute_gradphitilde();

  // then compute the  vtherm.(gradphi+gradphitilde) in fourier space
  compute_vtherm_dot_gradphi_plus_gradphitilde();
  

}

void FixGridModelH::post_force()
{

  double qx,qy,qz,q2;
  
  for (int i = 0; i < grid->chempot->Nz(); i++) {
    qz = grid->qzs[i];
    for (int j = 0; j < grid->chempot->Ny(); j++) {
      qy = grid->qys[j];
      for (int k = 0; k < grid->chempot->Nx(); k++) {



	qx = domain->dqx()*k;

	q2 = qx*qx + qy*qy + qz*qz;

	
	(*grid->chempot)(i,j,k) 
	  += temp/volFH*gamma*(*grid->laplace_phi)(i,j,k);
      }
    }
  }
}

void FixGridModelH::pre_final_integrate()
{


  fftw_execute(grid->forward_chempot);

  compute_vdet();

  didnotintegrate = true;
  return;
}


void FixGridModelH::final_integrate()
{

  conjugate->update();
  didnotintegrate = false;
  
  return;
}


void FixGridModelH::post_final_integrate()
{

  fftw_execute(grid->backward_phi);

  if (didnotintegrate) {
    double factor = grid->boxgrid[0]*grid->boxgrid[1]*grid->boxgrid[2];
    
    (*grid->phi) /= factor;
  }
  
  return;
}



void FixGridModelH::compute_vtherm_dot_gradphi()
{


  int localNx = vtherm_dot_gradphi->Nx();
  int localNy = vtherm_dot_gradphi->Ny();
  int localNz = vtherm_dot_gradphi->Nz();
  
  for (int i = 0; i < localNz; i++) 
    for (int j = 0; j < localNy; j++) 
      for (int k = 0; k < localNx; k++)
	(*vtherm_dot_gradphi)(i,j,k) = (*vtherm[0])(i,j,k)*(*gradphi[0])(i,j,k)
	  +(*vtherm[1])(i,j,k)*(*gradphi[1])(i,j,k)
	  +(*vtherm[2])(i,j,k)*(*gradphi[2])(i,j,k);


  fftw_execute(forward_vtherm_dot_gradphi);

}


void FixGridModelH::compute_vtherm_dot_gradphi_plus_gradphitilde()
{


  int localNx = vtherm_dot_gradphi->Nx();
  int localNy = vtherm_dot_gradphi->Ny();
  int localNz = vtherm_dot_gradphi->Nz();
  
  for (int i = 0; i < localNz; i++) 
    for (int j = 0; j < localNy; j++) 
      for (int k = 0; k < localNx; k++)
	(*vtherm_dot_gradphi)(i,j,k)
	  = (*vtherm[0])(i,j,k)*((*gradphi[0])(i,j,k)+(*gradphitilde[0])(i,j,k))
	  +(*vtherm[1])(i,j,k)*((*gradphi[1])(i,j,k)+(*gradphitilde[1])(i,j,k))
	  +(*vtherm[2])(i,j,k)*((*gradphi[2])(i,j,k)+(*gradphitilde[2])(i,j,k));


  fftw_execute(forward_vtherm_dot_gradphi);

}



void compute_gradphitilde()
{
  double normalization = 1.0/(grid->ft_boxgrid[0]*grid->ft_boxgrid[1]*grid->ft_boxgrid[2]);

  const int local0start = grid->ft_phi->get_local0start();
  const int globalNy = grid->ft_boxgrid[1];
  const int globalNz = grid->ft_boxgrid[2];
  
    
  std::complex<double> idqx(0,domain->dqx());
  std::complex<double> idqy(0,domain->dqy());
  std::complex<double> idqz(0,domain->dqz());
  
  double l,m,n;
  
  
  std::complex<double> tmp;
  
  for (int i = 0; i < grid->ft_phi->Nz(); i++) {
    
    if (i + local0start > globalNz/2) 
      l = -globalNz + i + local0start;
    else
      l = i + local0start;
    
    for (int j = 0; j < grid->ft_phi->Ny(); j++) {
      
      if (j > globalNy/2)
	m = -globalNy + j;
      else
	m = j;
      
      for (int k = 0; k < grid->ft_phi->Nx(); k++) {
	
	n = k;
	tmp = ft_noise(i,j,k)+(ft_phi(i,j,k)-dt*ft_vtherm_dot_gradphi(i,j,k))*normalization;
	
	(*grid->ft_gradphitilde[0])(i,j,k) = tmp*idqx*n;
	// SWAP HERE SINCE FFTW IS DOING A TRANSPOSED FT
	(*grid->ft_gradphitilde[1])(i,j,k) = tmp*idqz*l; 
	(*grid->ft_gradphitilde[2])(i,j,k) = tmp*idqy*m;
	
	
      }
    }
  }

  for (int i = 0; i < 3; i++)
    fftw_execute(backward_gradphitilde[i]);
  
}
