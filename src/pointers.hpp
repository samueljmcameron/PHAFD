#ifndef PHAFD_POINTERS_HPP
#define PHAFD_POINTERS_HPP


#include <mpi.h>
#include "phafd.hpp"

namespace PHAFD_NS {


class Pointers {
public:

  Pointers(PHAFD *ptr) : phafd(ptr), atoms(ptr->atoms),domain(ptr->domain),
			 grid(ptr->grid),commbrick(ptr->commbrick),
			 neighbor(ptr->neighbor),groups(ptr->groups),
			 input(ptr->input),atomfixes(ptr->atomfixes),
			 gridfixes(ptr->gridfixes),pairs(ptr->pairs),
			 world(ptr->world)  {};
  
protected:


  PHAFD *phafd;

  

  std::unique_ptr<Atom> &atoms;
  std::unique_ptr<Domain> &domain;
  std::unique_ptr<Grid> &grid;
  std::unique_ptr<CommBrick> &commbrick;
  std::unique_ptr<Neighbor> &neighbor;
  std::vector<std::unique_ptr<Group>> &groups;
  std::unique_ptr<Input> &input;
  std::vector<std::unique_ptr<FixAtom>> &atomfixes;
  std::vector<std::unique_ptr<FixGrid>> &gridfixes;
  std::vector<std::unique_ptr<Pair>> &pairs;

  MPI_Comm &world;
};
}


#endif
