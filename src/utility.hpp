#ifndef PHAFD_UTILITY_HPP
#define PHAFD_UTILITY_HPP

#include <cstdint>
#include <mpi.h>
#include <string>
#include <vector>
#include <map>
#include "phafd.hpp"
#include <complex>


namespace fftwArr {
  template<typename>
  class array3D;
}

namespace PHAFD_NS {

  namespace utility {
    std::vector<std::string> split_line(std::string&);
    void replacePercentages(std::string &, int);
    
    void convertVariables(std::string &,
			  std::map<std::string, std::string> const&);

  
    void check_MPI_duplicates(const std::vector<int> &,MPI_Comm ,int ,int ,
			      std::string);

    int make_unique_seed(int,const MPI_Comm &, int, int);

    int find_brackets(std::string &);

    int find_index(std::string id, const std::vector<std::string> &);

    template <class T>
    void type_of_output(std::string, std::string ,T *);
    
    void find_array_component(const std::string &,
			      PHAFD *,
			      fftwArr::array3D<double>* ,
			      std::string check_output="");

    void find_array_component(const std::string &,
			      PHAFD *,
			      fftwArr::array3D<std::complex<double>>*,
			      std::string check_output="");
    
  }

union ubuf {
  double d;
  int64_t i;
  ubuf(const double &arg) : d(arg) {}
  ubuf(const int64_t &arg) : i(arg) {}
  ubuf(const int &arg) : i(arg) {}
};
}
#endif
