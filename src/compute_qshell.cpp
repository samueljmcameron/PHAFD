#include "compute_qshell.hpp"

#include <cmath>
#include <iostream>
#include "utility.hpp"
#include "grid.hpp"
#include "domain.hpp"
#include "fix.hpp"
#include "comm_brick.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;

#define INARR(i,j,k) (input_array[(k)+((i)*localNy+(j))*localNx])


ComputeQshell::ComputeQshell(PHAFD *phafd) : Compute(phafd) {
  vector = true;
  gather_flag = false;
}

ComputeQshell::~ComputeQshell()
{
}

void ComputeQshell::init(const std::vector<std::string> &v_line,
			 bool add_to_names) {

  Compute::init(v_line,add_to_names);


  std::string arrname = v_line.at(1);

  double q = std::stod(v_line.at(2));
  double qwidth = std::stod(v_line.at(3));
  qlo = q-qwidth/2;
  qhi = q+qwidth/2;
  //throw std::runtime_error("Invalid argument in compute qshell.");

  if (v_line.size() == 5) {
    if (v_line.at(4) == "gather")
      gather_flag = true;
  } else if (v_line.size() > 5)
    throw std::runtime_error("Invalid argument in compute qshell");

  std::string prefix = arrname.substr(0,2);
  std::string id = arrname.substr(2);

  int input_nc;
  if (prefix == "c_") {

    int index = utility::find_index(id,Compute::NAMES);
    compute = computes.at(index).get();

    local0start = compute->local0start;
    localNz = compute->localNz;
    localNy = compute->localNy;
    localNx = compute->localNx;

    input_array = compute->array.data();

    input_nc = compute->numberofcomponents;
    
    
    
  } else if (prefix == "f_") {
    
    int index = utility::find_index(id,Fix::NAMES);
    fix = fixes.at(index).get(); 

    local0start = fix->local0start;
    localNz = fix->localNz;
    localNy = fix->localNy;
    localNx = fix->localNx;

    input_array = fix->array.data();

    input_nc = fix->numberofcomponents;
  } else
    throw std::runtime_error("Invalid input array in "
			     "compute/ft/spherical/bin.");

  if (input_nc > 1)
    // throws if component of array of length > 1 is not specified
    throw std::runtime_error("Cannot have multi-component array"
			     + id
			     + std::string(" in compute qshell."
					   ));

  
  // initialise list of indices 
  initialise();
  local_counts = indices.size()/3;
  local_array.resize(local_counts);

  MPI_Allreduce(&local_counts,&global_counts,1,
		MPI_INT,MPI_SUM,world);
  
  if (gather_flag) {
    global_counts_array.resize(commbrick->nprocs);
    
    MPI_Allgather(&local_counts,1,MPI_INT,
		  global_counts_array.data(),
		  1,MPI_INT,world);
    
    
    displacements.resize(commbrick->nprocs);
    
    displacements.at(0) = 0;
    for (int i = 1; i < commbrick->nprocs; i++)
      displacements.at(i)
	= global_counts_array.at(i-1) + displacements.at(i-1);
    array.resize(global_counts);  
  } else
    array.resize(local_counts);





  
}


void ComputeQshell::initialise()
{
  double dqx(domain->dqx());
  double dqy(domain->dqy());
  double dqz(domain->dqz());

  double l,m,n;

  double q;

  int ibin;


  int globalNy = grid->ft_boxgrid[1];
  int globalNz = grid->ft_boxgrid[2];


  for (int i = 0; i < localNz; i++) {
    
    if (i + local0start > globalNz/2) 
      l = -globalNz + i + local0start;
    else
      l = i + local0start;


    for (int j = 0; j < localNy; j++) {
      
      if (j > globalNy/2)
	m = -globalNy + j;
      else
	m = j;
      
      for (int k = 0; k < localNx; k++) {

	n = k;
	q = sqrt(dqx*n*dqx*n + dqz*m*dqz*m + dqy*l*dqy*l);
	if (q >= qlo && q < qhi) {
	  indices.push_back(i);
	  indices.push_back(j);
	  indices.push_back(k);
	}


      }

    }

  }

}


void ComputeQshell::start_of_step()
{
  if (this_step) {
    if (compute != nullptr)	
      compute->this_step = true;
    else if (fix != nullptr)
      fix->this_step = true;
  }
}

void ComputeQshell::in_fourier()
{
  if (!this_step) return;
  if (gather_flag)
    in_fourier_templated<true>();
  else
    in_fourier_templated<false>();
}
template < bool GATHER>
void ComputeQshell::in_fourier_templated()
{
  int count = 0;

  int i,j,k;

  if (GATHER) {
    for (auto & ar : local_array) {
      i = indices[count++];
      j = indices[count++];
      k = indices[count++];
      
      ar = INARR(i,j,k);
    }
    

      MPI_Allgatherv(local_array.data(),local_counts,MPI_DOUBLE,
		     array.data(),global_counts_array.data(),
		     displacements.data(),MPI_DOUBLE,world);

  } else
    for (auto & ar : array) {
      i = indices[count++];
      j = indices[count++];
      k = indices[count++];
      
      ar = INARR(i,j,k);
    }

}

void ComputeQshell::end_of_step()
{
  if (!this_step) return;
  if (compute != nullptr)	
    compute->this_step = false;
  else if (fix != nullptr)
    fix->this_step = false;
}
