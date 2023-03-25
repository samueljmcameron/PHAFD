#ifndef PHAFD_DUMP_HPP
#define PHAFD_DUMP_HPP

#include <set>
#include <string>
#include <fstream>
#include <Eigen/Core>

#include "pointers.hpp"

namespace psPDE {
  template<typename>
  class fftw_MPI_3Darray;
}

namespace PHAFD_NS {

class Dump : protected Pointers
{
public:
  Dump(PHAFD *);


  void init(const std::vector<std::string> &);

  void setup();
  void start_of_step();
  void require_calculations();
  void write_collection_header();
  void write_collection_footer();
  void write_collection_middle();
  
  inline static std::vector<std::string> NAMES;

  int every;

  std::string name;
  

protected:
  std::string base_name; // base name of file without extension, e.g. vtkfiles/field_p0
  std::string nopath_base_name; // name of file without path or extension, e.g. field_p0
  std::string collection_name; // collection name, .e.g vtkfiles/field_p0.pvd
  std::string instance_name; // base name plus time and extension, e.g. vtkfiles/field_p0_0.vti
  std::string nopath_instance_name; // above but without path, e.g. field_p0_0.vti

  std::string fext; // file extension of per time files (e.g. ".vti" for grid, ".vtp" for atom)
  
  std::string dump_type; // either atom, grid, or ftgrid


  std::set<std::string> attributes;

private:
  void create_instance_name();

  void write_atom_timestep();
  void write_grid_timestep();

  void process_attribute_name(std::fstream &, const std::string &);
  void write_ascii_data(std::fstream &, const std::string &,
			Eigen::Ref<Eigen::Matrix3Xd>);

  template <typename T>
  void write_ascii_data(std::fstream &, const std::string &,
			std::vector<T>,int);

  void append_binary_data(std::fstream &,
			  psPDE::fftw_MPI_3Darray<double> * ); 

  void append_binary_data(std::fstream &,const double *);
  
  std::string which_atoms;

  unsigned int bytelength;

  int precision;
  
};

}
#endif
