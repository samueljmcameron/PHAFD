#include "compute_ft_spherical_bin.hpp"

#include <cmath>
#include <iostream>
#include "utility.hpp"
#include "grid.hpp"
#include "fix.hpp"
#include "domain.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;


ComputeFtSphericalBin::ComputeFtSphericalBin(PHAFD *phafd) : Compute(phafd) {
  vector = true;
  numberofcomponents = 3;
}

ComputeFtSphericalBin::~ComputeFtSphericalBin()
{
}

void ComputeFtSphericalBin::init(const std::vector<std::string> &v_line,
				 bool add_to_names) {

  Compute::init(v_line,add_to_names);

  compute = nullptr;
  fix = nullptr;

  std::string arrname = v_line.at(1);

  
  try {
    input_component = utility::find_brackets(arrname);
  } catch (const std::runtime_error &e) {
    input_component = -1;
  }

  


  std::string prefix = arrname.substr(0,2);
  std::string id = arrname.substr(2);
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

  if (input_component == -1) {
    if (input_nc > 1)
      // throws if component of array of length > 1 is not specified
      throw std::runtime_error("Must specify component of multi-component "
			       "array " + id
			       + std::string(" in compute/ft/spherical/bin."
					     ));
    else
      input_component = 0;
  }

  nbins = std::stoi(v_line.at(2));
  double qmax = std::stod(v_line.at(3));
  dqbins = qmax/nbins;
  std::cout << "dqbins = " << dqbins << std::endl;
  output.resize(nbins);


  // array will have qs in 0..nbins-1, counts in nbins..2*nbins-1,
  //  and actual binned data in 2*nbins-1..3*nbins-1
  array.resize(nbins*numberofcomponents);

  counts.resize(nbins);

  global_counts.resize(nbins);


  // set counts and global counts to zero
  std::fill(counts.begin(),counts.end(),0);
  global_counts = counts;

  
  loop<1>();

  MPI_Allreduce(counts.data(),global_counts.data(),nbins,
		MPI_INT,MPI_SUM,world);


  for (int ibin = 0; ibin < nbins; ibin++)
    array.at(ibin) = dqbins*(0.5+ibin);
  for (int ibin = 0; ibin < nbins; ibin++)
    array.at(ibin+nbins) = global_counts[ibin];
  
  
}


void ComputeFtSphericalBin::start_of_step()
{
  if (this_step) {
    if (compute != nullptr)	
      compute->this_step = true;
    else if (fix != nullptr)
      fix->this_step = true;
  }
}
void ComputeFtSphericalBin::end_of_step()
{
  if (!this_step) return;

  std::fill(output.begin(),output.end(),0);

  loop<0>();
  MPI_Allreduce(output.data(),array.data()+2*nbins,nbins,
		MPI_DOUBLE,MPI_SUM,world);

  for (int ibin = 0; ibin < nbins; ibin++) {
    array[ibin+2*nbins] = array[ibin+2*nbins]/global_counts[ibin];
  }

  if (compute != nullptr)	
    compute->this_step = false;
  else if (fix != nullptr)
    fix->this_step = false;


}
template <int Tp_COUNT>
void ComputeFtSphericalBin::loop()
{


  
  double dqx(domain->dqx());
  double dqy(domain->dqy());
  double dqz(domain->dqz());

  int globalNy = grid->ft_boxgrid[1];
  int globalNz = grid->ft_boxgrid[2];


  double l,m,n;

  double q;

  int ibin;
  int index = input_component;
  
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

	ibin = static_cast<int>(q/dqbins);
	
	
	if (ibin < nbins) {
	  if (Tp_COUNT)
	    counts[ibin] += 1;
	  else
	    output[ibin] += input_array[index];
	}

	index += input_nc;

      }

    }

  }
  
  
  return;
  
}


