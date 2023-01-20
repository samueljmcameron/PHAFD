
#include "utility.hpp"
//#include <algorithm>
#include <exception>
#include <set>



void PHAFD::check_MPI_duplicates(const std::vector<int> & vec,
				 MPI_Comm comm,int id,int mpi_size,
				 std::string vecname)
/*
  check that atom IDs are not duplicated across processors, and that all vectors are filled.
  input vector should not have duplicate values (on same processor) but does not need to be
  sorted.
*/ 
{

  

  std::vector<int> number_to_check;

  if (id == 0)
    number_to_check.resize(mpi_size);

  
  int sendsize = vec.size();
  
  // gather all sendsizes to root
  MPI_Gather(&sendsize,1,MPI_INT,number_to_check.data(),1,MPI_INT,0,comm);
  

  
  std::vector<int> vec_to_check;
  std::vector<int> displs;


  
  if (id == 0) {
    
    displs.resize(mpi_size);
    
    int totalsize = 0;
    int count = 0;
    for (int num : number_to_check) {
      displs[count++] = totalsize;
      totalsize += num;
    }

    vec_to_check.resize(totalsize);
  }

  MPI_Gatherv(vec.data(),sendsize,MPI_INT,vec_to_check.data(),number_to_check.data(),
	      displs.data(),MPI_INT,0,comm);
  

  int intersections = 0;
  
  if (id == 0) {

    // check that there are no duplicates

    std::set<int> vec_set(vec_to_check.begin(),vec_to_check.end());

    intersections = vec_to_check.size()-vec_set.size();

    

  }
  
  MPI_Bcast(&intersections,1,MPI_INT,0,comm);
  
  if (intersections) throw std::runtime_error(std::to_string(intersections)
					      + std::string(" duplicate ") + vecname
					      + std::string(" across processors."));
  
  return;  
}
