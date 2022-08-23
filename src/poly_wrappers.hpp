#ifndef POLYWRAPPERS_HPP
#define POLYWRAPPERS_HPP

#include "beadrodpmer/single_tether.hpp"
#include "beadrodpmer/double_tether.hpp"
#include "beadrodpmer/no_tether.hpp"

#include <string>

class SingleTetherWrap : public BeadRodPmer::SingleTether {
public:

  int nuc_index;
  std::string name;
  std::vector<double> dFdX;

  SingleTetherWrap(const std::vector<std::string> & splitvec,
		   int nuc_index,std::string name)
    : BeadRodPmer::SingleTether(splitvec),nuc_index(nuc_index),
      name(name),dFdX{0,0,0} {};
  
  ~SingleTetherWrap() {};


};

class NoTetherWrap : public BeadRodPmer::NoTether {
public:

  int nuc_index;
  std::string name;
  std::vector<double> dFdX;
  
  NoTetherWrap(const std::vector<std::string> & splitvec,
	       int nuc_index,std::string name)
    : BeadRodPmer::NoTether(splitvec),nuc_index(nuc_index),
      name(name),dFdX{0,0,0} {};

  
  ~NoTetherWrap() {};


};

class DoubleTetherWrap : public BeadRodPmer::DoubleTether {
public:

  int nuc_index;
  std::string name;
  std::vector<double> dFdX;

  DoubleTetherWrap(const std::vector<std::string> & splitvec,
		   int nuc_index,std::string name)
    : BeadRodPmer::DoubleTether(splitvec),nuc_index(nuc_index),
      name(name),dFdX{0,0,0} {};

  
  ~DoubleTetherWrap() {};


};


#endif
