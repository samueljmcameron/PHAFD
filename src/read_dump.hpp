#ifndef PHAFD_READ_DUMP_HPP
#define PHAFD_READ_DUMP_HPP

#include <set>
#include <string>
#include <fstream>
#include <Eigen/Core>

#include "pointers.hpp"

namespace fftwArr {
  template<typename>
  class array3D;
}

namespace PHAFD_NS {

class ReadDump : protected Pointers
{
public:
  ReadDump(PHAFD *);
  ~ReadDump();

  void init(const std::vector<std::string> &);

  void process_attributes();

private:
  std::string filename; // filename to be read
  
  std::string dump_type; // either atom or grid

  std::set<std::string> attributes; // what properties are to be read in

  std::unique_ptr<fftwArr::array3D<double>> dummy_fftwArr;



  void process_atom_attributes();
  void process_grid_attributes();

  void read_ascii_data(const std::string &,
		       Eigen::Ref<Eigen::Matrix3Xd>);

  
  void read_ascii_data(const std::string &,std::vector<int> &);

  void read_binary_data(std::fstream &,
			fftwArr::array3D<double> * ); 
  
};

}
#endif
