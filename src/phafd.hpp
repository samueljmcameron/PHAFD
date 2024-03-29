#ifndef PHAFD_PHAFD_HPP
#define PHAFD_PHAFD_HPP


#include <vector>
#include <memory>
#include <mpi.h>
#include <fstream>
#include <map>

namespace PHAFD_NS {


class PHAFD 
{
  void create();
public:
  
  PHAFD(MPI_Comm,int, char **);
  ~PHAFD();


  /*
  class Atom * atoms;
  class Domain * domain;
  class Grid * grid;
  class CommBrick * commbrick;
    */
  std::unique_ptr<class Atom> atoms;
  std::unique_ptr<class Domain> domain;
  std::unique_ptr<class Grid> grid;
  std::unique_ptr<class CommBrick> commbrick;
  std::unique_ptr<class Neighbor> neighbor;
  std::unique_ptr<class Input> input;
  std::unique_ptr<class Integrate> integrate;
  std::vector<std::unique_ptr<class Group>> groups;
  std::vector<std::unique_ptr<class Fix>> fixes;
  std::vector<std::unique_ptr<class Pair>> pairs;
  std::vector<std::unique_ptr<class Compute>> computes;
  std::vector<std::unique_ptr<class Dump>> dumps;
  
  
  MPI_Comm world;

  std::unique_ptr<std::fstream> inputfile;
  
  std::unique_ptr<std::map<std::string,std::string>> varmap;
  
};
}
#endif
