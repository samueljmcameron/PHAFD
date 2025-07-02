#ifndef PHAFD_COMPUTE_HPP
#define PHAFD_COMPUTE_HPP


#include <set>
#include <complex>

#include "pointers.hpp"

namespace fftwArr {
  template<typename>
  class array3D;
}


namespace PHAFD_NS {

class Compute : protected Pointers
{
public:
  Compute(PHAFD *);

  std::string name;
  std::vector<double> array;
  virtual void init(const std::vector<std::string> &);
  virtual void in_fourier() = 0;
  virtual void end_of_step() = 0;

  void start_of_step();
  
  inline static std::vector<std::string> NAMES;
  bool per_grid;
  bool per_ftgrid;
  bool per_atom;
  bool scalar;
  bool vector;

  bool clusterscomputed;

  bool this_step;
  std::set<std::string> dump_callers;


  std::vector<fftwArr::array3D<double> *> realFFTWarray;
  std::vector<fftwArr::array3D<std::complex<double>> *> complexFFTWarray;

  
  //int localNx,localNy,localNz;

  int numberofcomponents; // number of components in the array (e.g. 1 for scalar, 3 for vector, etc.)
  
};

}
#endif
