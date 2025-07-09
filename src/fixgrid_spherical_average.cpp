
#include "grid.hpp"
#include "domain.hpp"
#include "fixgrid_spherical_average.hpp"
#include "fftw_arr/array3d.hpp"



FixGridSphericalAverage::FixGridSphericalAverage(PHAFD *phafd)
  : Fix(*phafd),vector(true);
{

}

void FixGridSphericalAverage::init(const std::vector<std::string> &v_line)
{

  Fix::init(v_line);


  nbins = std::stoi(v_line.at(1));
  dqbins = std::stod(v_line.at(2));

  output.resize(nbins);
  array.resize(nbins);

  counts.resize(nbins);

  global_counts.resize(nbins);
  global_counts = counts;
  
  fftwArr::array3D<std::complex<double>> *temporary_array;
  
  if(grid->ft_phi != nullptr) 
    temporary_array = grid->ft_phi.get();
  else if(grid->ft_velocity[0] != nullptr)
    temporary_array = grid->ft_velocity[0].get();
  else
    throw std::runtime_error("Need ft_phi or ft_velocity allocated!");

  local0start = temporary_array->get_local0start();
  localNz = temporary_array->Nz();
  localNy = temporary_array->Ny();
  localNx = temporary_array->Nx();


  globalNy = grid->ft_boxgrid[1];
  globalNz = grid->ft_boxgrid[2];

}


void FixGridSphericalAverage::setup()
{

  // set counts and global counts to zero
  std::fill(counts.begin(),counts.end(),0);
  global_counts = counts;
  
  loop<1>(nullptr);



  MPI_Allreduce(counts.data(),global_counts.data(),nbins,
		MPI_DOUBLE,MPI_SUM,world);
} 

void FixGridSphericalAverage::end_of_step()
{

  std::fill(output.begin(),output.end(),0);


  loop<0>(ft_phi.get());

  MPI_Allreduce(output.data(),array.data(),nbins,
		MPI_DOUBLE,MPI_SUM,world);

  for (int ibin = 0; ibin < nbins; ibin++) {
    array[ibin] = array[ibin]/global_counts[ibin];
  }
  
}

template <int Tp_COUNT>
void FixGridSphericalAverage::loop(const fftwArr::array3D<double> *fftw_arr)
{


  
  double dqx(domain->dqx());
  double dqy(domain->dqy());
  double dqz(domain->dqz());
  
  double l,m,n;

  double q;

  int ibin;
  
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


	q = sqrt(dqx*n*dqx*n + dqz*m*dqz*m + dqy*l*dqy*l);

	ibin = q/dqbin;


	if (ibin < nbins) {
	  if (Tp_COUNT)
	    counts[ibin] += 1;
	  else:
	    output[ibin] += (*fftw_arr)(i,j,k);
	}

      }

    }

  }
  
  
  return;
  
}

