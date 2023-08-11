#include "compute_grid_clusters.hpp"

#include <cmath>

#include "grid.hpp"
#include "comm_brick.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;

#define INTARR(i,j,k) (intarray[(k)+((i)*Ny+(j))*Nx])
#define RIGHTS(j,k) (rights[(j)*Ny+(k)])
#define LEFTS(j,k) (lefts[(j)*Ny+(k)])
#define MIN(A,B) ((A) < (B) ? (A) : (B))


ComputeGridClusters::ComputeGridClusters(PHAFD *phafd) : Compute(phafd) {
  per_grid = true;
  clusterscomputed = true;
}


void ComputeGridClusters::init(const std::vector<std::string> &v_line) {

  Compute::init(v_line);
  
  
  std::string arrname = v_line.at(1);
  int globalNz;

  // iterate over different arrays that are available eventually... for now just assume it's phi

  if (arrname != "phi")
    throw std::runtime_error("Must use phi as array for clustering.");
  
  fftw3_arr = grid->phi.get();

  if (fftw3_arr == nullptr)
    throw std::runtime_error("Cannot compute complex of array " + arrname
			     + std::string(", not allocated."));


  threshold = std::stod(v_line.at(2));

  condition = v_line.at(3);

  if (condition != "lt" && condition != "gt")
    throw std::runtime_error("Must specify lt or gt for compute grid cluster.");
  

  Nz = fftw3_arr->Nz();
  Ny = fftw3_arr->Ny();
  Nx = fftw3_arr->Nx();
  
  array.resize(Nx*Ny*(Nz+2));

  intarray.resize(Nx*Ny*Nz);
  
  lefts.resize(Nx*Ny);
  rights.resize(Nx*Ny);


  
}

void ComputeGridClusters::end_of_step()
{

  if (!this_step) return;

  int *clusint;

  int local_start = fftw3_arr->get_local0start();

  if (condition == "lt") {
  
    for (int i = 0; i < Nz; i++) {
      for (int j = 0; j < Ny; j++) {
	for (int k = 0; k < Nx; k++) {
	  if ((*fftw3_arr)(i,j,k) < threshold)
	    INTARR(i,j,k) = i+local_start + (j*Nx + k)*Ny;
	  else
	    INTARR(i,j,k) = -1;
	}
      }
    }
  } else {

    for (int i = 0; i < Nz; i++) {
      for (int j = 0; j < Ny; j++) {
	for (int k = 0; k < Nx; k++) {
	  if ((*fftw3_arr)(i,j,k) > threshold)
	    INTARR(i,j,k) = i+local_start + (j*Nx + k)*Ny;
	  else
	    INTARR(i,j,k) = -1;
	}
      }
    }
  }
    
  // send bulk grid data to the neighboring border data
  send_to_neighbors();


  int change,anychange;
  bool done;

  while (true) {

    // send neighboring border data back to the bulk process


    send_to_neighbors();//    receive_from_neighbors();

    change = 0;


    MPI_Barrier(world);
    
    while (true) {
      done = true;

      for (int i = 0; i < Nz; i++) {
	for (int j = 0; j < Ny; j++) {
	  for (int k = 0; k < Nx; k++) {

	    
	    // if not part of cluster, don't care
	    if (INTARR(i,j,k) == -1) continue;
	    
	    // check over nearest neighbors

	    for (int iindex = i-1; iindex <= i+1; iindex += 2) {
	      if (iindex == -1)
		clusint = &LEFTS(j,k);
	      else if (iindex == Nz)
		clusint = &RIGHTS(j,k);
	      else
		clusint = &INTARR(iindex,j,k);

	      
	      if (*clusint == -1 || *clusint == INTARR(i,j,k)) continue;

	      
	      INTARR(i,j,k) = *clusint = MIN(INTARR(i,j,k),*clusint);

	      
	      done = false;
	    }

	    
	    

	    for (int jindex = j-1; jindex <= j+1; jindex += 2) {
	      if (jindex == -1)
		clusint = &INTARR(i,Ny-1,k);
	      else if (jindex == Ny)
		clusint = &INTARR(i,0,k);
	      else
		clusint = &INTARR(i,jindex,k);

	      if (*clusint == -1 || *clusint == INTARR(i,j,k)) continue;
	      
	      INTARR(i,j,k) = *clusint = MIN(INTARR(i,j,k),*clusint);
	      
	      done = false;
	    }


	    for (int kindex = k-1; kindex <= k+1; kindex += 2) {
	      if (kindex == -1)
		clusint = &INTARR(i,j,Nx-1);
	      else if (kindex == Nx)
		clusint = &INTARR(i,j,0);
	      else
		clusint = &INTARR(i,j,kindex);

		  
	      if (*clusint == -1 || *clusint == INTARR(i,j,k)) continue;

	      
	      INTARR(i,j,k) = *clusint = MIN(INTARR(i,j,k),*clusint);
	      
	      
	      done = false;

	    }

	    


	  }
	}
      }

      if (!done) change = 1;
      if (done) break;

    }
    
    MPI_Allreduce(&change,&anychange,1,MPI_INT,MPI_MAX,world);
    if (!anychange) break;

  }


  int count = 0;
  if (commbrick->me != 0) {
    for (int j = 0; j < Ny; j++)
      for (int k = 0; k < Nx; k++)
	array[count++] = LEFTS(j,k);
  }
  
  for (int i = 0; i < Nz; i++) {
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nx; k++) {
	array[count++] = INTARR(i,j,k);
      }
    }
  }

  if (commbrick->me != commbrick->nprocs-1) {
    for (int j = 0; j < Ny; j++)
      for (int k = 0; k < Nx; k++)
	array[count++] = RIGHTS(j,k);
  }
  
  
}

void ComputeGridClusters::send_to_neighbors()
{
  int nprocs = commbrick->nprocs;
  int me = commbrick->me;
  int recvid,sendid;
  
  // send to the right/recv from the left
  if (me == 0) {
    recvid = nprocs-1;
    sendid = me+1;

  } else if (me == nprocs-1) {
    recvid = me-1;
    sendid = 0;
    
  } else {
    recvid = me -1;
    sendid = me+1;
  }

  MPI_Sendrecv(&INTARR(Nz-1,0,0),Nx*Ny,MPI_INT,sendid,0,
	       lefts.data(),Nx*Ny,MPI_INT,recvid,0,world,MPI_STATUS_IGNORE);
  
  
  // now send to left/recv from right
  if (me == 0) {
    recvid = me+1;
    sendid = nprocs-1;

  } else if (me == nprocs-1) {
    recvid = 0;
    sendid = me-1;
    
  } else {
    recvid = me+1;
    sendid = me -1;
  }

  MPI_Sendrecv(&INTARR(0,0,0),Nx*Ny,MPI_INT,sendid,0,
	       rights.data(),Nx*Ny,MPI_INT,recvid,0,world,MPI_STATUS_IGNORE);
  
}
/*
// must receive in order somehow..
void ComputeGridClusters::receive_from_neighbors()
{
  int nprocs = commbrick->nprocs;
  int me = commbrick->me;
  int recvid,sendid;
  
  // send to the right/recv from the left
  if (me == 0) {
    recvid = nprocs-1;
    sendid = me+1;

  } else if (me == nprocs-1) {
    recvid = me-1;
    sendid = 0;
    
  } else {
    recvid = me -1;
    sendid = me+1;
  }



  for (int j = 0; j < Ny; j++) {
    for (int k = 0; k < Nx; k++) {

      RIGHTS(j,k) = INTARR(Nz-1,j,k);
      LEFTS(j,k) = INTARR(0,j,k);

    }
  }
  
  MPI_Sendrecv(rights.data(),Nx*Ny,MPI_INT,sendid,0,
	       &INTARR(0,0,0),Nx*Ny,
	       MPI_INT,recvid,0,world,MPI_STATUS_IGNORE);  


  
  // now send to left/recv from right
  if (me == 0) {
    recvid = me+1;
    sendid = nprocs-1;

  } else if (me == nprocs-1) {
    recvid = 0;
    sendid = me-1;
    
  } else {
    recvid = me+1;
    sendid = me -1;
  }



  MPI_Sendrecv(lefts.data(),Nx*Ny,MPI_INT,sendid,0,
	       &INTARR(Nz-1,0,0),Nx*Ny,
	       MPI_INT,recvid,0,world,MPI_STATUS_IGNORE);
  
  
}
*/
