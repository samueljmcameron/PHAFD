
#ifndef PHAFD_FIXATOM_SEMIFLEXIBLE_HPP
#define PHAFD_FIXATOM_SEMIFLEXIBLE_HPP


#include "fixatom.hpp"

#include <memory>

namespace PHAFD_NS {

template <typename T>
class FixAtomSemiFlexible : public FixAtom {
public:
  FixAtomSemiFlexible(PHAFD *);

  virtual void init(const std::vector<std::string> &) override;
  
  virtual void setup() override;

  virtual void initial_integrate() override;
  virtual void final_integrate() override;
private:

  std::vector<std::unique_ptr<T>> pmers;
  std::vector<int> nbeads;

  void check_bonds();
  
};

}

#endif
