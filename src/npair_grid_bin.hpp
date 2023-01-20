

#ifndef PHAFD_NPAIR_GRID_BIN_H
#define PHAFD_NPAIR_GRID_BIN_H

#include "npair.hpp"

namespace PHAFD {

  
class NPairGridBin : public NPair {
 public:
  NPairGridBin();
  void build(class NeighList *, const Atom &) override;
  void build(class NeighList *, const Atom &,
	     const psPDE::Grid &) override;
  
};

}    // namespace LAMMPS_NS

#endif
