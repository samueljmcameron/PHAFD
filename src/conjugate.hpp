#ifndef PHAFD_CONJUGATE_HPP
#define PHAFD_CONJUGATE_HPP

#include <string>
#include <complex>
#include <vector>

#include "domain.hpp"
#include "grid.hpp"
#include "smatrix.hpp"
#include "pointers.hpp"

namespace fftwArr {
  template <typename T>
  class array3D;
}

namespace PHAFD_NS {
class Conjugate : protected Pointers
{
  
public:
  Conjugate(PHAFD *);


  virtual void reset_dt(double);
  void setup();
  void update();

  virtual void readCoeffs(const std::vector<std::string> &) = 0;


private:
  
  const int nblocks;
  int me,nprocs;

protected:
  
  fftwArr::array3D<std::complex<double>> *ft_array;

  
  double complexprefactor, realprefactor, sqrtdt;
  

  std::vector<double> qys,qzs;

  
  double dt;
  
  int localNy,localNz;

  // arrays for transferring via mpi
  std::vector<sMatrix<std::complex<double>>> blocks;
  std::vector<sMatrix<std::complex<double>>> lines;

  std::vector<MPI_Request> block_requests;
  std::vector<MPI_Request> line_requests;


  
  

  bool equal_flag; // true if all processors have same number of grid points
  int divider;     // processor whose domain has middle line (at globalNz/2) in it


  virtual void complex_update(int,int,int) = 0;
  virtual void real_update(int,int,int) = 0;
  virtual void origin_update() = 0;

  void set_qs();
  
  void single(int);
  
  void first(int);
  void lefthalf(int);
  void righthalf(int);
  void middle_odd(int);
  void middle_even(int);
  void last(int);

  void line_single(int);
  
  void line_first(int);
  void line_lefthalf(int);
  void line_righthalf(int);
  void line_middle_odd(int);
  void line_middle_even(int);
  void line_last(int);



  
  void onehalf_equal(int, std::string);
  void onehalf_unequal(int,std::string);
  

  void middle_odd_equal(int);
  void middle_odd_unequal(int);
  void middle_even_equal(int);
  void middle_even_unequal(int);


  void line_onehalf_equal(int,std::string);
  void line_onehalf_unequal(int,std::string);
  
  void line_lefthalf_unequal(int);
  void line_righthalf_unequal(int);
  void line_middle_odd_equal(int);
  void line_middle_odd_unequal(int);
  void line_middle_even_equal(int);
  void line_middle_even_unequal(int);

  

  
  
};
};

#endif
