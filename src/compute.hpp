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

  virtual void init(const std::vector<std::string> &,bool);
  virtual void in_fourier() = 0;
  virtual void start_of_step();
  virtual void end_of_step() = 0;

  
  inline static std::vector<std::string> NAMES;
  bool per_grid; // this is true when the compute can output onto a grid dump
  //                    i.e. data stored in realFFTWarray of size (Nx,Ny,Nz)
  bool per_ftgrid; // this is true when the compute can output onto a ftgrid dump
  //                    i.e. data stored in realFFTWarray but of size (Nx/2+1,Ny,Nz)
  bool per_atom; // this is true when the compute can output to a (lammps) dump
  //                    i.e. data stored in array<double> with numberofcomponents*Natoms
  bool scalar;   // this is true when there is single scalar output
  bool vector;   // this is true when there is  array with numberofcomponents*M but M
  //                    M is not Natoms (so can't output to dump) - example is outputting
  //                    spherical binning in Fourier space

  bool clusterscomputed;

  bool this_step; // usually this is set externally by other fixes/dumps
  std::set<std::string> dump_callers;


  /*
    vectors here are used to pass info to another compute or fix only,
    they are not output to e.g. dump files
  */
  std::vector<fftwArr::array3D<double> *> realFFTWarray;
  std::vector<fftwArr::array3D<std::complex<double>> *> complexFFTWarray;

  // array is what is output to dump files.
  std::vector<double> array;
  int numberofcomponents; // number of components in the array (e.g. 1 for scalar, 3 for vector, etc.)
  int localNx,localNy,localNz,local0start;  
};

}
#endif
