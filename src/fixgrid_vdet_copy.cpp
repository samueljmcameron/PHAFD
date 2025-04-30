#include <algorithm>

#include "utility.hpp"


#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "integrate.hpp"
#include "fixgrid_vdet.hpp"
#include "conjugate_volfrac.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;


FixGridVdet::FixGridVdet(PHAFD *phafd) : Fix(phafd) {};



void FixGridVdet::init(const std::vector<std::string> &v_line)
/*
  v_line should take form:
  fixname,seedx,seedy,seedz,viscosity,temperature
 */
{

  Fix::init(v_line);

  std::vector<std::string> new_v_line;

  

  std::array<std::string,3> xys_list = {"x","y","z"};
  std::array<int,3> seeds;

  // build string vector for concentration


  for (int iarg = 0; iarg < 3; iarg++)  {
    new_v_line.clear();
    new_v_line.push_back("v"+std::string(xyz_list.at(iarg)));
    seeds.at(iarg) = std::stoi(v_line.at(iarg+1));
    seeds.at(iarg) = utility::make_unique_seed(seeds.at(iarg),
					       world,commbrick->me,commbrick->nprocs);
    new_v_line.push_back(seeds.at(iarg));
    new_v_line.push_back(v_line.at(4));
    new_v_line.push_back(v_line.at(5));
    conjugate_noise.at(iarg) = std::make_unique<ConjugateNoise>(phafd);
    conjugate_noise.at(iarg)->readCoeffs(new_v_line);
  }


  qys = conjugate_vnoise.at(0)->qys.data();
  qzs = conjugate_vnoise.at(0)->qzs.data();

  viscosity = conjugate_vnoise.at(0)->damping;
  
  ft_Znoise_x = ft_Znoise[0].data();
  ft_Znoise_y = ft_Znoise[1].data();
  ft_Znoise_z = ft_Znoise[2].data();


  
}

  
void FixGridVdet::setup()
{

}

void FixGridVdet::reset_dt()
{

  dt = integrate->dt;

  for (int i = 0; i < 3; i++)
    conjugate_vnoise.at(i)->reset_dt(dt);
  
}


void FixGridVdet::pre_final_integrate()
{

  compute_chempot_grad_phi();
  
}



void FixGridVdet::compute_chempot_gradphi() {

  
  for (int i = 0; i < grid->chempot->Nz(); i++) 
    for (int j = 0; j < grid->chempot->Ny(); j++) 
      for (int k = 0; k < grid->chempot->Nx(); k++) {
	(*velocity[0])(i,j,k) = (*chempot)(i,j,k)*(*gradphi[0])(i,j,k);
	(*velocity[1])(i,j,k) = (*chempot)(i,j,k)*(*gradphi[1])(i,j,k);
	(*velocity[2])(i,j,k) = (*chempot)(i,j,k)*(*gradphi[2])(i,j,k);
      }


  


	
}


// setting the values of vdet in fourier space, transposes are accounted for here.

void FixGridVdet::set_vdet(int i, int j, int k) {

  double qx,qy,qz,q2,Txx,Txy,Txz,Tyy,Tyz,Tzz;


  qz = qzs[i];
  qy = qys[j];
  qx = domain->dqx()*k;
  
  q2 = qx*qx + qy*qy + qz*qz;

  if (q2 == 0) {
    (*ft_vdet[0])(i,j,k) = 0.0;
    (*ft_vdet[1])(i,j,k) = 0.0;
    (*ft_vdet[2])(i,j,k) = 0.0;
  } else {



    Txx = (1.0-qx*qx/q2)/(q2*viscosity);
    Txy = (-qx*qy/q2)/(q2*viscosity);  
    Txz = (-qx*qz/q2)/(q2*viscosity);
    Tyy = (1-qy*qy/q2)/(q2*viscosity);
    Tyz = (-qy*qz/q2)/(q2*viscosity);
    Tzz = (1-qz*qz/q2)/(q2*viscosity);


    // need to do a swap here (qy <-> qz) since computing transposed
    // fourier functions, the below LOOKS LIKE IT HAS BUGS BUT IT DOES NOT!!

    (*ft_vdet[0])(i,j,k) = (Txx*(*ft_Znoise_x)(i,j,k)+Txz*(*ft_Znoise_y)(i,j,k)
			      + Txy*(*ft_Znoise_z)(i,j,k))/dt;
    (*ft_vdet[1])(i,j,k) = (Txz*(*ft_Znoise_x)(i,j,k)+Tzz*(*ft_Znoise_y)(i,j,k)
			      + Tyz*(*ft_Znoise_z)(i,j,k))/dt;
    (*ft_vdet[2])(i,j,k) = (Txy*(*ft_Znoise_x)(i,j,k)+Tyz*(*ft_Znoise_y)(i,j,k)
			      + Tyy*(*ft_Znoise_z)(i,j,k))/dt;    

  }

  
  return ;

}
