#ifndef PHAFD_INPUT_HPP
#define PHAFD_INPUT_HPP

#include "pointers.hpp"

namespace PHAFD_NS {

class Input : protected Pointers {
public:
  Input(PHAFD *);

  ~Input();

  void read();


};
  
}
#endif
