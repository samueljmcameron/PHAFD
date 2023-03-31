#ifndef PHAFD_READ_ATOMS_HPP
#define PHAFD_READ_ATOMS_HPP

#include <string>
#include <fstream>


#include "pointers.hpp"

namespace PHAFD_NS {


class ReadAtoms : protected Pointers {
public:

  static inline int SUCCESS = 0;
  static inline int FORMAT_ERROR = 1;
  static inline int NOFILE = 2;

  ReadAtoms(PHAFD *);

  int read_file(const std::string &);


private:



  std::string fname;
  std::ifstream datafile;

  
  int create_atoms(int);

  

};


}

#endif
