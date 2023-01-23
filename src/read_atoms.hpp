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




  std::ifstream datafile;


  std::vector<std::string> styles; // names of the particles styles (size 1 if not hybrid)
  std::vector<int> particlesPerStyle; // # of particles of each style
  std::vector<int> global_particlesPerStyle;

  
  std::vector<int> atomsPerStyle; // # of atoms of each style
  std::vector<int> global_atomsPerStyle;

  std::vector<int> atomsPerPolymer;
  std::vector<int> global_atomsPerPolymer;


  std::vector<int> global_particlesPerProc;

  
  


  int read_perStyle();


  int reset_styles();
  int reset_particlesPerStyle();
  int reset_atomsPerStyle();
  
  int create_atoms();

  

};


}

#endif
