#ifndef PHAFD_COMPUTE_PAIR_HPP
#define PHAFD_COMPUTE_PAIR_HPP


#include "compute.hpp"

namespace PHAFD_NS {

class ComputePair : public Compute
{
public:
  ComputePair(PHAFD *);


  virtual void init(const std::vector<std::string> &) override;
  virtual void in_fourier() override {};
  virtual void end_of_step() override;

  
private:

  class Pair *pair;
};
  

}
#endif
