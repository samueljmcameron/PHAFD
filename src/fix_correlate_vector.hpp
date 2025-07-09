
#ifndef PHAFD_CORRELATE_VECTOR_HPP
#define PHAFD_CORRELATE_VECTOR_HPP


#include "fix.hpp"

#include <memory>

namespace PHAFD_NS {

class Compute;

class FixCorrelateVector : public Fix {
public:
  FixCorrelateVector(PHAFD *);

  virtual void init(const std::vector<std::string> &,
		    bool) override;
  
  virtual void setup() override;

  virtual void initial_integrate() override {};
  virtual void post_force() override {};
  virtual void pre_final_integrate() override {};
  virtual void final_integrate() override {};
  virtual void post_final_integrate(bool /* flag */) override {};
  

  virtual void start_of_step() override;
  virtual void end_of_step() override;

private:

  int every,repeat,freq;
  int Nstores;

  int elements_in_storage,circular_index;

  double *input_array;
  int input_nc, input_array_size;

  std::string output_type;
  
  Compute *compute;
  Fix *fix;

  
  std::vector<int> Nsamples;

  std::vector<double> storage_array;

  double *get_storage_section(int);
  void add_to_array(int,double *, double *);
  
};

}

#endif
