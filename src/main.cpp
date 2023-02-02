#include <mpi.h>
#include <fftw3-mpi.h>
#include "phafd.hpp"
#include "input.hpp"
#include "domain.hpp"
#include "utility.hpp"

int main(int argc, char **argv)
{
  // imagine that the system is split into mpi_size rectangles

  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Info info = MPI_INFO_NULL;

  int ierr = MPI_Init(NULL,NULL);

  fftw_mpi_init();

  {

    std::string line = "boxdims 100.0 100.0 100.0 boxorigin -50.0 -50.0 -50.0";
    auto v_line = PHAFD_NS::utility::split_line(line);

    
    auto phafd = new PHAFD_NS::PHAFD(comm,argc,argv);
    phafd->input->read();

    
    delete phafd;


    
  }
  fftw_mpi_cleanup();

  
  ierr = MPI_Finalize();
  return 0;

}

