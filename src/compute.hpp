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

  bool this_step;
  std::set<std::string> dump_callers;


  /*
    realFFTWarray is a vector of arrays which can be used in two ways:
    (1) to store an array for future fix/compute calculations
    (2) to store an array which is output to a dump file, EITHER agrid
        or ftgrid type!

    In either case, the array size is EITHER has indices up to (Nx,Ny,Nz)
    or  (Nx/2+1,Ny,Nz) (in theory it could be used to store any sized
    3D array, but these are two typical cases). This might be confusing as
    most other arrays of this type are always (Nx,Ny,Nz). But this is necessary
    for outputting ftgrid type data to dump, since this ftgrid data must be
    of the latter size.
  */
  std::vector<fftwArr::array3D<double> *> realFFTWarray;


  /*
    complexFFTWarray is a vector of arrays which is just used to store
    complex data types, but it CANNOT be output to dump files as only
    real arrays can be output to dump files (by design).
  */
  std::vector<fftwArr::array3D<std::complex<double>> *> complexFFTWarray;

  
  //int localNx,localNy,localNz;

  int numberofcomponents; // number of components in the array (e.g. 1 for scalar, 3 for vector, etc.)
  
};

}
#endif
