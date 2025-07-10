#ifndef PHAFD_READ_VTP_HPP
#define PHAFD_READ_VTP_HPP

#include <vector>
#include <string>
#include <fstream>
#include <Eigen/Core>

#include "pointers.hpp"

namespace fftwArr {
  template<typename>
  class array3D;
}

namespace PHAFD_NS {

class ReadVTP : protected Pointers
{
public:
  ReadVTP(PHAFD *);
  ~ReadVTP();

  void init(const std::vector<std::string> &);

  void read();

private:
  std::string filename; // filename to be read
  std::string directory; // filename to be read

  std::vector<std::string> v_line_for_read_dump;
  
  std::vector<std::string> list_of_files;
  std::ifstream vtpfile, pvtifile;
  std::string dump_type; // either atom or grid

  std::vector<std::string> attributes; // what properties are to be read in

  bool no_padding;
  void read_preamble(const std::string &,const std::string &,
		     std::ifstream &);
  std::string read_vtp_for_pvti_name();
  void create_file_list();
  std::string read_in_vti();
  std::string::size_type check_pvd_extension();
  int64_t check_instance_name(const std::string &);
  int64_t check_pinstance_name(const std::string &);

  int64_t Nfirst, Nlast, every, skip;
  
  
};

}
#endif
