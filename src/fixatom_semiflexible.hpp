
#ifndef PHAFD_FIXATOM_SEMIFLEXIBLE_HPP
#define PHAFD_FIXATOM_SEMIFLEXIBLE_HPP


#include "fixatom.hpp"

#include <memory>

namespace PHAFD {

template <typename T>
class FixAtomSemiFlexible : public FixAtom {
public:
  FixAtomSemiFlexible(const Group &, std::string,Atom &,Domain &);

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
