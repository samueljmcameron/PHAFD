#ifndef POLYWRAPPERS_HPP
#define POLYWRAPPERS_HPP

#include "beadrodpmer/single_tether.hpp"
#include "beadrodpmer/double_tether.hpp"
#include "beadrodpmer/no_tether.hpp"

#include <string>

class SingleTetherWrap : public BeadRodPmer::SingleTether {
public:

  int number,nuc_index;
  std::string name;
  std::vector<std::vector<double>> dFdX_is;

  SingleTetherWrap(const std::vector<std::string> & splitvec,
		   int number,int nuc_index,std::string name)
    : BeadRodPmer::SingleTether(splitvec),number(number),nuc_index(nuc_index),
      name(name) {};
  
  ~SingleTetherWrap() {};


};

class NoTetherWrap : public BeadRodPmer::NoTether {
public:

  int number,nuc_index;
  std::string name;
  std::vector<std::vector<double>> dFdX_is;
  
  NoTetherWrap(const std::vector<std::string> & splitvec,
	       int number,int nuc_index,std::string name)
    : BeadRodPmer::NoTether(splitvec),number(number),nuc_index(nuc_index),
      name(name) {};

  
  ~NoTetherWrap() {};


};

class DoubleTetherWrap : public BeadRodPmer::DoubleTether {
public:

  int number,nuc_index;
  std::string name;
  std::vector<std::vector<double>> dFdX_is;

  DoubleTetherWrap(const std::vector<std::string> & splitvec,
		   int number,int nuc_index,std::string name)
    : BeadRodPmer::DoubleTether(splitvec),number(number),nuc_index(nuc_index),
      name(name) {};

  
  ~DoubleTetherWrap() {};


};


#endif
