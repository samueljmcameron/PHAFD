

#ifndef PHAFD_NPAIR_GRID_BIN_H
#define PHAFD_NPAIR_GRID_BIN_H

#include "npair.hpp"

namespace PHAFD_NS {

  
class NPairGridBin : public NPair {
 public:
  NPairGridBin(PHAFD *);

  void build(class NeighList *) override;
  
};

}    // namespace LAMMPS_NS

#endif
