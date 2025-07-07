
#ifndef PHAFD_OUTPUT_FILE_HPP
#define PHAFD_OUTPUT_FILE_HPP


#include "fix.hpp"

#include <memory>

namespace PHAFD_NS {

class Compute;

class FixOutputFile : public Fix {
public:
  FixOutputFile(PHAFD *);

  virtual void init(const std::vector<std::string> &) override;
  
  virtual void setup() override;

  virtual void initial_integrate() override {};
  virtual void post_force() override {};
  virtual void pre_final_integrate() override {};
  virtual void final_integrate() override {};
  virtual void post_final_integrate(bool /* flag */) override {};
  

  virtual void start_of_step() override;
  virtual void end_of_step() override;

private:

  int every;

  double *input_array;
  int input_nc, input_array_size;

  std::string filename,output_type;
  
  Compute *compute;
  Fix *fix;

  std::vector<double> output;
};

}

#endif
