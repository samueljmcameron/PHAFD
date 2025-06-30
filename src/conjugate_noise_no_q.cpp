
#include "conjugate_noise_no_q.hpp"
#include "grid.hpp"
#include <cmath>
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;

ConjugateNoiseNoQ::ConjugateNoiseNoQ(PHAFD *phafd)
  : ConjugateNoise(phafd)
{}

void ConjugateNoiseNoQ::complex_update(int i , int j, int k)
{

  (*ft_array)(i,j,k).real(complexprefactor*sqrtdt*real_dist(gen));
  (*ft_array)(i,j,k).imag(complexprefactor*sqrtdt*real_dist(gen));

  
  return;
  
}


void ConjugateNoiseNoQ::real_update(int i, int j, int k)
{


  (*ft_array)(i,j,k).real(realprefactor*sqrtdt*real_dist(gen));
  
  return;
}
