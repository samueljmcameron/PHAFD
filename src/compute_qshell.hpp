#ifndef PHAFD_COMPUTE_QSHELL_HPP
#define PHAFD_COMPUTE_QSHELL_HPP

#include <vector>
#include "compute.hpp"


namespace fftwArr {
  template<typename>
  class array3D;
}

namespace PHAFD_NS {

class Fix;
  
class ComputeQshell : public Compute
{
public:
  ComputeQshell(PHAFD *);
  ~ComputeQshell();


  virtual void init(const std::vector<std::string> &,bool) override;
  virtual void in_fourier() override;
  virtual void start_of_step() override;
  virtual void end_of_step() override;

  
  
private:

  void initialise();
  template <int Tp_COUNT> void loop();
  template < bool GATHER> void in_fourier_templated();
  std::vector<int> indices;


  bool gather_flag;
  Compute *compute;
  Fix *fix;

  double *input_array;

  std::vector<double> local_array;
  std::vector<int> global_counts_array;
  std::vector<int> displacements;

  double qlo, qhi;


  int local_counts, global_counts;

};
  

}
#endif
