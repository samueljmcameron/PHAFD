#ifndef PHAFD_COMM_BRICK_HPP
#define PHAFD_COMM_BRICK_HPP

#include <vector>
#include <array>

#include "pointers.hpp"

/*==================================================================
  Structure expected is a 3D system where the x and y only have one
  processor, but z is split into nproc processors
  ==================================================================*/

namespace PHAFD_NS {
class CommBrick : protected Pointers {
public:

  CommBrick(PHAFD *);
  ~CommBrick();
  void borders();
  void forward_comm();
  void reverse_comm();
  
  void setup(double);

  std::vector<double> buf_send,buf_recv; // flat vector of size 3*numatoms to send/recv
  int me,nprocs;
  
private:

  double zprd,zinterval,zlo,cutghost;
  const int xyswaps;

  std::vector<double> subloz,subhiz;
  int nswap;

  int size_forward,size_reverse,size_border;
  
  std::array<double,3> maxneed;
  
  std::vector<bool> pbc_flag;
  std::vector<std::array<double,3>> pbc;
  std::vector<double> slablo, slabhi;

  std::vector<int> sendproc,recvproc,firstrecv,sendnum,recvnum;
  std::vector<int> size_forward_recv,size_reverse_recv,size_reverse_send;
  std::vector<std::vector<int>> sendlists;
  std::vector<int> labels;
  std::vector<std::vector<double>> pbcz;

};
}
#endif
