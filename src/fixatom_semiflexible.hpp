
#ifndef PHAFD_FIXATOM_SEMIFLEXIBLE_HPP
#define PHAFD_FIXATOM_SEMIFLEXIBLE_HPP


#include "fix.hpp"

#include <memory>

namespace PHAFD_NS {

template <typename T>
class FixAtomSemiFlexible : public Fix {
public:
  FixAtomSemiFlexible(PHAFD *);

  virtual void init(const std::vector<std::string> &) override;
  
  virtual void setup() override;
  virtual void initial_integrate() override;
  virtual void post_force() override {};
  virtual void pre_final_integrate() override {};
  virtual void final_integrate() override;
  virtual void post_final_integrate() override {};

  virtual void end_of_step() override {};

private:

  std::vector<std::unique_ptr<T>> pmers;
  std::vector<int> nbeads;

  void check_bonds();

  bool collect_this_step;
  
};

}

#endif
