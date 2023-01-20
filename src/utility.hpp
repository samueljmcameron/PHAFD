#ifndef PHAFD_UTILITY_HPP
#define PHAFD_UTILITY_HPP

#include <cstdint>
#include <mpi.h>
#include <string>
#include <vector>

namespace PHAFD {

  void check_MPI_duplicates(const std::vector<int> &,MPI_Comm ,int ,int ,
  			    std::string);
union ubuf {
  double d;
  int64_t i;
  ubuf(const double &arg) : d(arg) {}
  ubuf(const int64_t &arg) : i(arg) {}
  ubuf(const int &arg) : i(arg) {}
};
}
#endif
