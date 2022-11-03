#ifndef POLYMERWRAP_HPP
#define POLYMERWRAP_HPP

#include "beadrodpmer/polymer.hpp"
#include "beadrodpmer/initialise.hpp"

#include <string>
#include <vector>
#include <memory>

class PolymerWrap 
{ 
private:
  int read_index;
public:

  int number,nuc_index;
  std::string name;
  std::vector<std::vector<double>> dFdX_is;

  bool read;

  PolymerWrap(std::unique_ptr<BeadRodPmer::Polymer> _poly,
		   int number,int nuc_index,std::string _name)
    : polymer(std::move(_poly)),number(number),nuc_index(nuc_index),read(false)
  {
    std::string::size_type start = _name.find("read_");

    if (start != std::string::npos) {
      _name.erase(start,5);
      read = true;
      std::string::size_type last = _name.find_last_not_of("0123456789");
      read_index = std::stoi(_name.substr(last+1));
      _name.erase(last+1);
    }

    name = _name + std::to_string(number);
    
  };

  int get_read_index() {
    if (!read)
      throw std::runtime_error("trying to access read index for polymer which doesn't read in data.");
    return read_index;
  };


  void setup()
  {
    polymer->compute_tangents_and_friction();
    polymer->set_Hhat();
    polymer->set_dCdlambda();
    polymer->set_G();
  }

  int single_step(double t, double dt,int itermax, int numtries,
		  bool throw_exception)
  {
    return polymer->single_step(t,dt,dFdX_is,itermax,numtries,throw_exception);
  }

  void initialise_atoms(const std::vector<std::string> & splitvec)
  {
    BeadRodPmer::Initialise::init_atoms(splitvec,
					polymer->atoms,polymer->initspringK,
					polymer->initdt,polymer->inittolerance,
					polymer->equilibration_steps);
  }

  const std::vector<BeadRodPmer::Atom> & get_atoms() const
  { return polymer->atoms;};

  const std::vector<int> & get_nuc_beads() const 
  { return polymer->nuc_beads;};

  const BeadRodPmer::Polymer & get_polymer() const { return *polymer;};
  BeadRodPmer::Polymer & get_polymer() { return *polymer;};
private:
  std::unique_ptr<BeadRodPmer::Polymer> polymer;

};

#endif
