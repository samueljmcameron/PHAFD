#include "comm_brick.hpp"
#include "conjugate.hpp"
#include "grid.hpp"
#include "domain.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;
Conjugate::Conjugate(PHAFD *phafd)
  : Pointers(phafd),nblocks(4),ft_array(nullptr)
{
  me = commbrick->me;
  nprocs = commbrick->nprocs;

}


void Conjugate::reset_dt(double timestep) {

  dt = timestep;
  return;
}


void Conjugate::setup()
{

  localNy = ft_array->Ny();
  localNz = ft_array->Nz();
  
  
  
  std::vector<int> localNzs(nprocs);
  
  for (int i = 0; i < nprocs; i++) {
    
    if (i == me) {
      localNzs[i] = localNz;
      
    }
    
    MPI_Bcast(&localNzs[i],1,MPI_INT,i,world);

  }
  

  set_qs();

  
  // I'm assuming that the complex array blocks will always be split in a certain
  // configuration (see logic in if statement on next line).
  for (int i = 1; i < nprocs; i++ ) {
    if (localNzs.at(nprocs-i) == 0 ) {
      std::cout << "EXITING!!! number of processors chosen such that me = " << nprocs-i
		<< " has array size zero." << std::endl;
      MPI_Abort(world,1);
    } else if (localNzs.at(0) < localNzs.at(i)) {
      std::cout << "EXITING!!! Block size of process (me = " << i
		<< ") is larger than block size of initial process me = 0." << std::endl;
      MPI_Abort(world,1);
    } else if (localNzs.at(i-1) != localNzs.at(i) && i != nprocs-1) {
      std::cout << "EXITING!!! Block size of process (me = " << i
		<< ") is not the same size as the previous block." << std::endl;
      MPI_Abort(world,1);
    }
    
  }


  // now set up matrices which store data to communicate via MPI
  int x = localNzs.at(0);
  int y = localNzs.at(nprocs-1);

  
  block_requests.resize(nblocks);
  line_requests.resize(nblocks);


  if (x == y) {

    equal_flag = true;
    divider = nprocs/2;
    
    blocks.push_back(sMatrix<std::complex<double>>(1,localNy/2));
    blocks.push_back(sMatrix<std::complex<double>>(x-1,localNy/2));
    blocks.push_back(sMatrix<std::complex<double>>(1,localNy/2));
    blocks.push_back(sMatrix<std::complex<double>>(x-1,localNy/2));

    lines.push_back(sMatrix<std::complex<double>>(1,1));
    lines.push_back(sMatrix<std::complex<double>>(x-1,1));
    lines.push_back(sMatrix<std::complex<double>>(1,1));
    lines.push_back(sMatrix<std::complex<double>>(x-1,1));
    
  } else {
    
    equal_flag = false;
    divider = (nprocs-1)/2;

    blocks.resize(nblocks);
    lines.resize(nblocks);
    
    if (me == 0 || me == nprocs-1) {
      blocks[0] = sMatrix<std::complex<double>>(y,localNy/2);
      blocks[2] = sMatrix<std::complex<double>>(y,localNy/2);

      lines[0] = sMatrix<std::complex<double>>(y,1);
      lines[2] = sMatrix<std::complex<double>>(y,1);

      
    } else {
      blocks[0] = sMatrix<std::complex<double>>(y+1,localNy/2);
      blocks[2] = sMatrix<std::complex<double>>(y+1,localNy/2);

      lines[0] = sMatrix<std::complex<double>>(y+1,1);
      lines[2] = sMatrix<std::complex<double>>(y+1,1);
      
    }

    blocks[1] = sMatrix<std::complex<double>>(x-y-1,localNy/2);
    blocks[3] = sMatrix<std::complex<double>>(x-y-1,localNy/2);

    lines[1] = sMatrix<std::complex<double>>(x-y-1,1);
    lines[3] = sMatrix<std::complex<double>>(x-y-1,1);
    
  }    

}


/* Set the values of arrays for qx, qy, and qz */
void Conjugate::set_qs()
{


  const int local0start = ft_array->get_local0start();

  const int globalNy = grid->ft_boxgrid[1];
  const int globalNz = grid->ft_boxgrid[2];

  const double dy = domain->dqy();
  const double dz = domain->dqz();
  
  qys.resize(localNy);
  qzs.resize(localNz);


  for (int i = 0; i < localNz; i++)
    if (i + local0start > globalNz/2) 
      qzs[i] = -dz*(globalNz - i -local0start);
    else 
      qzs[i] = dz*(i+local0start);


  for (int j = 0; j < localNy; j++)
    if (j > globalNy/2) 
      qys[j] = -dy * (globalNy - j);
    else
      qys[j] = dy * j;

  return;

}

void Conjugate::update()
{

  const int localNx = ft_array->Nx();  
  const int local0start = ft_array->get_local0start();
  
  const int globalNx = grid->ft_boxgrid[0];
  const int globalNy = grid->ft_boxgrid[1];
  const int globalNz = grid->ft_boxgrid[2];
  
  
  
  // update the bulk of the system which has no constraints
  
  for (int nz = 0; nz < localNz; nz ++) {
    
    for (int ny = 0; ny < localNy; ny ++) {
      
      for (int nx = 1; nx < localNx-1; nx++) {
	  complex_update(nz,ny,nx);
      }
    }
  }
  
  
  
  // update the two planes at nx = 0 and nx = Nx-1
  // which are constrained to ensure that phi remains real
    
  for (int nx = 0; nx < localNx; nx += localNx-1) {
    
    if (nprocs == 1) {
      
      single(nx);
      line_single(nx);
      
    } else {
      
      if (equal_flag && me == 0) {
	first(nx);
	line_first(nx);
      } else if (me < divider) {
	lefthalf(nx);
	line_lefthalf(nx);
      } else if (me == divider) {
	if (nprocs % 2 == 0) {
	  middle_even(nx);
	  line_middle_even(nx);
	} else {
	  middle_odd(nx);
	  line_middle_odd(nx);
	  }
      } else if (!equal_flag && me == nprocs-1) {
	last(nx);
	line_last(nx);
      } else {
	righthalf(nx);
	line_righthalf(nx);
      }
    }
  }
  
  
  // update the eight points where ft_noise must be real
  if (me == 0) {
    int nx = localNx-1;
    real_update(0,localNy/2,0);
    real_update(0,localNy/2,nx);
    real_update(0,0,nx);
    
    origin_update();
    
  }
  if (local0start <= globalNz/2
      && local0start + localNz -1 >= globalNz/2) {
    int nz = globalNz/2 - local0start;
    int nx = localNx-1;
    real_update(nz,0,0);
    real_update(nz,localNy/2,0);
    real_update(nz,localNy/2,nx);
    real_update(nz,0,nx);
  }
  
  return;
}







void Conjugate::single(int nx)
/* Only to be used when there is only one processor. */
{

  
  
  // write line at nz = 0  
  int nz = 0;
  for (int j = 1; j < localNy/2; j++) {
    
    complex_update(nz,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
    (*ft_array)(nz,localNy-j,nx) = std::conj((*ft_array)(nz,j,nx));
					     
  }
  
  
  // write bottom left square to top right square where possible
  for (int i = 1; i < localNz/2; i++) {
    
    for (int j = 1; j < localNy/2; j++) {
      
      complex_update(i,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      (*ft_array)(localNz-i,localNy-j,nx) = std::conj((*ft_array)(i,j,nx));
      
    }
  }
  // write bottom right square to top left square where possible
  for (int i = localNz/2+1; i < localNz; i++) {
    for (int j = 1; j < localNy/2; j++ ) {
      complex_update(i,j,nx);
      (*ft_array)(localNz-i,localNy-j,nx) = std::conj((*ft_array)(i,j,nx));
    }
  }
  
  
  // don't forget Nz/2 term!
  nz = localNz/2;
  
  
  for (int j = 1; j < localNy/2; j++) {
    
    complex_update(nz,j,nx);//1.0*(nz+complex_local) + 1i*(1.0*j);
    (*ft_array)(nz,localNy-j,nx) = std::conj((*ft_array)(nz,j,nx));
    
  }
  
  return;
}

void Conjugate::line_single(int nx)
/* Only to be used when there is only one processor. */
{


  
  int ny = 0;
  // write bottom left square to bottom right square where possible
  for (int i = 1; i < localNz/2; i++) {
    complex_update(i,ny,nx);//1.0*(i+complex_local)+1i*(1.0*j);
    (*ft_array)(localNz-i,ny,nx) = std::conj((*ft_array)(i,ny,nx));
      
  }
  ny = localNy/2;
  
  // write top right square to top left square where possible
  for (int i = localNz/2+1; i < localNz; i++) {
    complex_update(i,ny,nx);
    (*ft_array)(localNz-i,ny,nx) = std::conj((*ft_array)(i,ny,nx));
  }

  return;
  
}

  


void Conjugate::first(int nx)
/* only to be used when all processors have same
   number of data points on them. */
{
  


  // write the line at nz = 0;

  int nz = 0;
  for (int j = 1; j < localNy/2; j++) {
    
    complex_update(nz,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
    (*ft_array)(nz,localNy-j,nx) = std::conj((*ft_array)(nz,j,nx));
					     
  }
  
  // send from far left corner of bottom left square

  for (int i = 1; i < localNz; i++) {

    for (int j = 1; j < localNy/2; j++) {

      complex_update(i,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      //(ft_noise(i,j,nx) - mobility*q2*ft_nonlinear(i,j,nx)*dt + noise*dt )/(1+mobility*gamma*q2*q2);
      
      
      blocks[1](blocks[1].axis_size(0)-i,blocks[1].axis_size(1)-j) = std::conj((*ft_array)(i,j,nx));
      
    }
  }
  
  MPI_Isend(blocks[1].data(),blocks[1].size()*2,MPI_DOUBLE,nprocs-1,1,
	    world,&block_requests[0]);
  
  // receive from bottom right square 
  MPI_Recv(blocks[3].data(), blocks[3].size()*2,MPI_DOUBLE,nprocs-1,3,
	   world,MPI_STATUS_IGNORE);
  
  for (int i = 1; i < localNz; i++) {
    for (int j = localNy/2+1; j < localNy; j++) {
      
      (*ft_array)(i,j,nx) = blocks[3](i-1,j-blocks[3].axis_size(1));
      
      
    }
  }
  
  
  MPI_Wait(&block_requests[0],MPI_STATUS_IGNORE);
  
  return;

}

void Conjugate::last(int nx)
/* only to be used when the last processor has
   less data points than the rest. */

{
  


  
  // send from far right corner of bottom right square
  for (int i = 0; i < blocks[2].axis_size(0); i++) {

    for (int j = 1; j < localNy/2; j++) {
      
      
      complex_update(i,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      //(ft_noise(i,j,nx) - mobility*q2*ft_nonlinear(i,j,nx)*dt + noise*dt )/(1+mobility*gamma*q2*q2);
      
      
      blocks[2](blocks[2].axis_size(0)-i-1,blocks[2].axis_size(1)-j) = std::conj((*ft_array)(i,j,nx));
      
    }
  }
  
  MPI_Isend(blocks[2].data(),blocks[2].size()*2,MPI_DOUBLE,nprocs-me-1,2,
	    world,&block_requests[2]);
  
  // receive from bottom left square 1 <= y < Nz/2, z <= 1 < Ny/2
  MPI_Recv(blocks[0].data(), blocks[0].size()*2,MPI_DOUBLE,nprocs-me-1,0,
	   world,MPI_STATUS_IGNORE);
  
  for (int i = 0; i < blocks[0].axis_size(0); i++) {
    for (int j = localNy/2+1; j < localNy; j++) {
      
      (*ft_array)(i,j,nx) = blocks[0](i,j-blocks[0].axis_size(1));
      
      
    }
  }
  
  
  MPI_Wait(&block_requests[2],MPI_STATUS_IGNORE);
  
  return;

}


void Conjugate::lefthalf(int nx)
{
  if (equal_flag)
    onehalf_equal(nx,"left");
  else
    onehalf_unequal(nx,"left");
  return;
}

void Conjugate::righthalf(int nx)
{
  if (equal_flag)
    onehalf_equal(nx,"right");
  else
    onehalf_unequal(nx,"right");
  return;
}


void Conjugate::middle_odd(int nx)
{
  if (equal_flag)
    middle_odd_equal(nx);
  else
    middle_odd_unequal(nx);
  return;
}

void Conjugate::middle_even(int nx)
{
  if (equal_flag)
    middle_even_equal(nx);
  else
    middle_even_unequal(nx);
  return;
}


void Conjugate::line_first(int nx)
/* only to be used when all processors have same
   number of data points on them. */

{
  
  int ny = 0;
  
  // send from far left corner of bottom left square
  for (int i = 1; i < localNz; i++) {

    
    complex_update(i,ny,nx);//1.0*(i+complex_local)+1i*(1.0*j);
    //(ft_noise(i,j,nx) - mobility*q2*ft_nonlinear(i,j,nx)*dt + noise*dt )/(1+mobility*gamma*q2*q2);
      
      
    lines[1](lines[1].axis_size(0)-i,0) = std::conj((*ft_array)(i,ny,nx));
      
    
  }
  
  MPI_Isend(lines[1].data(),lines[1].size()*2,MPI_DOUBLE,nprocs-1,1,
	    world,&line_requests[0]);
  
  // receive from bottom right square 
  MPI_Recv(lines[3].data(), lines[3].size()*2,MPI_DOUBLE,nprocs-1,3,
	   world,MPI_STATUS_IGNORE);
  
  ny = localNy/2;
  for (int i = 1; i < localNz; i++) {
    (*ft_array)(i,ny,nx) = lines[3](i-1,0);
  }
  
  
  MPI_Wait(&line_requests[0],MPI_STATUS_IGNORE);
  
  return;

}


void Conjugate::line_last(int nx)
/* only to be used when the last processor has
   less data points than the rest. */
  
{
  
  int ny = localNy/2;
  
  // send from far right corner of bottom right square
  for (int i = 0; i < lines[2].axis_size(0); i++) {
    complex_update(i,ny,nx);//1.0*(i+complex_local)+1i*(1.0*j);
    //(ft_noise(i,j,nx) - mobility*q2*ft_nonlinear(i,j,nx)*dt + noise*dt )/(1+mobility*gamma*q2*q2);

    lines[2](lines[2].axis_size(0)-i-1,0) = std::conj((*ft_array)(i,ny,nx));

  }
  
  MPI_Isend(lines[2].data(),lines[2].size()*2,MPI_DOUBLE,nprocs-me-1,2,
	    world,&line_requests[2]);
  
  // receive from bottom left square 1 <= y < Nz/2, z <= 1 < Ny/2
  MPI_Recv(lines[0].data(), lines[0].size()*2,MPI_DOUBLE,nprocs-me-1,0,
	   world,MPI_STATUS_IGNORE);

  ny = 0;
  for (int i = 0; i < lines[0].axis_size(0); i++) {
    
    (*ft_array)(i,ny,nx) = lines[0](i,0);

  }
  
  MPI_Wait(&line_requests[2],MPI_STATUS_IGNORE);
  
  return;

}


void Conjugate::line_lefthalf(int nx)
{
  if (equal_flag)
    line_onehalf_equal(nx,"left");
  else
    line_onehalf_unequal(nx,"left");
  return;
}

void Conjugate::line_righthalf(int nx)
{
  if (equal_flag)
    line_onehalf_equal(nx,"right");
  else
    line_onehalf_unequal(nx,"right");
  return;
}


void Conjugate::line_middle_odd(int nx)
{
  if (equal_flag)
    line_middle_odd_equal(nx);
  else
    line_middle_odd_unequal(nx);
  return;
}

void Conjugate::line_middle_even(int nx)
{
  if (equal_flag)
    line_middle_even_equal(nx);
  else
    line_middle_even_unequal(nx);
  return;
}

void Conjugate::onehalf_equal(int nx,
			      std::string half)
/* only to be used when all processors have same
   number of data points on them. */
{

  int isend_0, jsend_1;
  int irec_0, jrec_1;

  if (half == "left") {
    isend_0 = 0;
    jsend_1 = 1;
    irec_0 = 2;
    jrec_1 = 3;
  } else if (half == "right") {
    isend_0 = 2;
    jsend_1 = 3;
    irec_0 = 0;
    jrec_1 = 1;
  } else { // should never get here
    std::cout << "fucked up in coding." << std::endl;
    return;
  }

  
  // single line at start of bottom left plane square

  int nz = 0;
  for (int j = 1; j < localNy/2; j++) {
    
    complex_update(nz,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
    //(ft_noise(i,j,nx) - mobility*q2*ft_nonlinear(i,j,nx)*dt + noise*dt )/(1+mobility*gamma*q2*q2);
    
    blocks[isend_0](nz,blocks[isend_0].axis_size(1)-j) = std::conj((*ft_array)(nz,j,nx));
    
  }
  

  
  MPI_Isend(blocks[isend_0].data(),blocks[isend_0].size()*2,MPI_DOUBLE,nprocs-me,isend_0,
	    world,&block_requests[isend_0]);
  
  for (int i = 1; i < localNz ; i++) {
    
    for (int j = 1; j < localNy/2; j++ ) {
      
      complex_update(i,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      blocks[jsend_1](blocks[jsend_1].axis_size(0)-i,blocks[jsend_1].axis_size(1)-j)
	= std::conj((*ft_array)(i,j,nx));
      
    }
  }
  
  MPI_Isend(blocks[jsend_1].data(),blocks[jsend_1].size()*2,MPI_DOUBLE,nprocs-me-1,jsend_1,
	    world,&block_requests[jsend_1]);
  
  
  
  // receive from bottom right square to fill in top left square
  // (1 <= y < Nz/2, Ny/2 < z < Ny)
  
  
  MPI_Irecv(blocks[irec_0].data(), blocks[irec_0].size()*2,MPI_DOUBLE,nprocs-me,irec_0,
	    world,&block_requests[irec_0]);
  
  MPI_Irecv(blocks[jrec_1].data(), blocks[jrec_1].size()*2,MPI_DOUBLE,nprocs-me-1,jrec_1,
	    world,&block_requests[jrec_1]);
  
  
  // now fill phi with blocks[jrec_1] and block4 received data in an efficient manner
  
  int b2flag=0;
  int b3flag = 0;
  
  
  // next block of code ensures that whichever block data is received first
  // is filled in first
  while (!(b2flag || b3flag)) {
    MPI_Test(&block_requests[irec_0],&b2flag,MPI_STATUS_IGNORE);
    MPI_Test(&block_requests[jrec_1],&b3flag,MPI_STATUS_IGNORE);
  }
  
  
  if (b2flag) { // if blocks[irec_0] is received first

    for (int j = localNy/2+1; j < localNy; j++) {
	
	(*ft_array)(nz,j,nx) = blocks[irec_0](nz,j-blocks[irec_0].axis_size(1));
	
    }
    
    MPI_Wait(&block_requests[jrec_1],MPI_STATUS_IGNORE);
    
    for (int i = 1; i < localNz; i++) {
      for (int j = localNy/2+1; j < localNy; j++) {
	
	(*ft_array)(i,j,nx) = blocks[jrec_1](i-1,j-blocks[jrec_1].axis_size(1));
	
      }
    }
  } else if (b3flag) { // if blocks[jrec_1] is received first

    for (int i = 1; i < localNz; i++) {
      for (int j = localNy/2+1; j < localNy; j++) {
	
	(*ft_array)(i,j,nx) = blocks[jrec_1](i-1,j-blocks[jrec_1].axis_size(1));
	
      }
    }
    
    MPI_Wait(&block_requests[irec_0],MPI_STATUS_IGNORE);
    
    
    
    for (int j = localNy/2+1; j < localNy; j++) {
      
      (*ft_array)(nz,j,nx) = blocks[irec_0](nz,j-blocks[irec_0].axis_size(1));
      
    }


  }
  
  MPI_Wait(&block_requests[isend_0],MPI_STATUS_IGNORE);    
  MPI_Wait(&block_requests[jsend_1],MPI_STATUS_IGNORE);    
  return;
}

void Conjugate::middle_odd_equal(int nx)
/* only to be used when all processors have same
   number of data points on them. */
{



  // send one to next processor

  int nz = 0;
  for (int j = 1; j < localNy/2; j++) {
    
    complex_update(nz,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
    //(ft_noise(i,j,nx) - mobility*q2*ft_nonlinear(i,j,nx)*dt + noise*dt )/(1+mobility*gamma*q2*q2);
    
    blocks[0](nz,blocks[0].axis_size(1)-j) = std::conj((*ft_array)(nz,j,nx));
    
  }
  

  
  MPI_Isend(blocks[0].data(),blocks[0].size()*2,MPI_DOUBLE,nprocs-me,0,
	    world,&block_requests[0]);


  // receive one from next processor

  MPI_Recv(blocks[2].data(), blocks[2].size()*2,MPI_DOUBLE,nprocs-me,2,
	   world,MPI_STATUS_IGNORE);

  for (int j = localNy/2+1; j < localNy; j++) {
    
    (*ft_array)(nz,j,nx) = blocks[2](nz,j-blocks[2].axis_size(1));
    
  }


  
  // write bottom left square to top right square where possible
  for (int i = 1; i < (blocks[1].axis_size(0)-1)/2+1; i++) {

    for (int j = 1; j < localNy/2; j++) {
      
      
      complex_update(i,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      (*ft_array)(localNz-i,localNy-j,nx) = std::conj((*ft_array)(i,j,nx));
      
    }
  }
  // write bottom right square to top left square where possible
  for (int i = (blocks[1].axis_size(0)-1)/2 + 2; i < localNz; i++) {
    for (int j = 1; j < localNy/2; j++ ) {
      complex_update(i,j,nx);
      (*ft_array)(localNz-i,localNy-j,nx) = std::conj((*ft_array)(i,j,nx));
    }
  }

  
  // don't forget Nz/2 term!
  nz = (blocks[1].axis_size(0)-1)/2+1;
  
  
  for (int j = 1; j < localNy/2; j++) {
    
    complex_update(nz,j,nx);//1.0*(nz+complex_local) + 1i*(1.0*j);
    (*ft_array)(nz,localNy-j,nx) = std::conj((*ft_array)(nz,j,nx));
    
  }

  MPI_Wait(&block_requests[0],MPI_STATUS_IGNORE);

  return;
  
}


void Conjugate::middle_even_equal(int nx)
/* only to be used when all processors have same
   number of data points on them. */
{




  // don't forget Nz/2 term!
  int nz = 0;
  
  
  for (int j = 1; j < localNy/2; j++) {
    
    complex_update(nz,j,nx);//1.0*(nz+complex_local) + 1i*(1.0*j);
    (*ft_array)(nz,localNy-j,nx) = std::conj((*ft_array)(nz,j,nx));
    
  }

  // fill and send bottom right square  
    
  for (int i = 1; i < localNz; i++) {
    for (int j = 1; j < localNy/2; j++) {
      
      complex_update(i,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      
      blocks[3](blocks[3].axis_size(0)-i,blocks[3].axis_size(1)-j)
	= std::conj((*ft_array)(i,j,nx));
    }
    
  }
  
  
  MPI_Isend(blocks[3].data(),blocks[3].size()*2,MPI_DOUBLE,nprocs-me-1,3,
	    world,&block_requests[3]);
  
  
  // receive from bottom left square, then fill in top right square
  
  MPI_Recv(blocks[1].data(), blocks[1].size()*2,MPI_DOUBLE,nprocs-me-1,1,
	    world,MPI_STATUS_IGNORE);
  

  for (int i = 1; i < localNz; i++) {
    for (int j = localNy/2+1; j < localNy; j++) {
      
      (*ft_array)(i,j,nx) = blocks[1](i-1,j-blocks[1].axis_size(1));
      
    }
  }

  MPI_Wait(&block_requests[3],MPI_STATUS_IGNORE);
  
  return;
  
}





void Conjugate::onehalf_unequal(int nx,
				std::string half)
/* only to be used when the last processor has
   less data points than the rest. */
{

  int irec_0,jrec_1,isend_0,jsend_1;

  if (half == "left") {
    isend_0 = 0;
    jsend_1 = 1;
    irec_0 = 2;
    jrec_1 = 3;
  } else if (half == "right") {
    isend_0 = 2;
    jsend_1 = 3;
    irec_0 = 0;
    jrec_1 = 1;
  } else {
    std::cout << "fucked up somehwere" << std::endl;
    return;
  }
    

  int loopstart = 0;
  if (me == 0) {
    loopstart = 1;
    // write the line at nz = 0;
    
    int nz = 0;
    for (int j = 1; j < localNy/2; j++) {
      
      complex_update(nz,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      (*ft_array)(nz,localNy-j,nx) = std::conj((*ft_array)(nz,j,nx));
      
    }
  }

  for (int i = loopstart; i < blocks[isend_0].axis_size(0)+loopstart; i++) {
    for (int j = 1; j < localNy/2; j++) {
      complex_update(i,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      blocks[isend_0](blocks[isend_0].axis_size(0)-i-(1-loopstart),blocks[isend_0].axis_size(1)-j) = std::conj((*ft_array)(i,j,nx));
    }
  }

  MPI_Isend(blocks[isend_0].data(),blocks[isend_0].size()*2,MPI_DOUBLE,nprocs-me-1,isend_0,
	    world,&block_requests[isend_0]);

  for (int i = blocks[isend_0].axis_size(0)+loopstart; i < localNz; i++) {
    for (int j = 1; j < localNy/2; j++) {
      complex_update(i,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      blocks[jsend_1](blocks[jsend_1].axis_size(0)-(i-blocks[isend_0].axis_size(0)+1-loopstart),blocks[jsend_1].axis_size(1)-j)
	= std::conj((*ft_array)(i,j,nx));
    }
  }
  
  MPI_Isend(blocks[jsend_1].data(),blocks[jsend_1].size()*2,MPI_DOUBLE,nprocs-me-2,jsend_1,
	    world,&block_requests[jsend_1]);
  MPI_Irecv(blocks[irec_0].data(), blocks[irec_0].size()*2,MPI_DOUBLE,nprocs-me-1,irec_0,
	    world,&block_requests[irec_0]);
  
  MPI_Irecv(blocks[jrec_1].data(), blocks[jrec_1].size()*2,MPI_DOUBLE,nprocs-me-2,jrec_1,
	    world,&block_requests[jrec_1]);

  int b2flag=0;
  int b3flag = 0;

  while (!(b2flag || b3flag)) {
    MPI_Test(&block_requests[irec_0],&b2flag,MPI_STATUS_IGNORE);
    MPI_Test(&block_requests[jrec_1],&b3flag,MPI_STATUS_IGNORE);
  }
  
  
  if (b2flag) { 
    for (int i = loopstart; i < blocks[irec_0].axis_size(0)+loopstart; i++) {
      for (int j = localNy/2+1; j < localNy; j++) {
	(*ft_array)(i,j,nx) = blocks[irec_0](i-loopstart,j-blocks[irec_0].axis_size(1));
      }
    }
    
    MPI_Wait(&block_requests[jrec_1],MPI_STATUS_IGNORE);
    
    for (int i = blocks[irec_0].axis_size(0)+loopstart; i < localNz; i++) {
      for (int j = localNy/2+1; j < localNy; j++) {
	
	(*ft_array)(i,j,nx) = blocks[jrec_1](i-(blocks[irec_0].axis_size(0)+loopstart),j-blocks[jrec_1].axis_size(1));
	
      }
    }
  } else if (b3flag) { // if blocks[jrec_1] is received first

    for (int i = blocks[irec_0].axis_size(0)+loopstart; i < localNz; i++) {
      for (int j = localNy/2+1; j < localNy; j++) {
	
	(*ft_array)(i,j,nx) = blocks[jrec_1](i-(blocks[irec_0].axis_size(0)+loopstart),j-blocks[jrec_1].axis_size(1));
	
      }
    }
    
    MPI_Wait(&block_requests[irec_0],MPI_STATUS_IGNORE);
    
    
    for (int i = loopstart; i < blocks[irec_0].axis_size(0)+loopstart; i++) {
      for (int j = localNy/2+1; j < localNy; j++) {
	
	(*ft_array)(i,j,nx) = blocks[irec_0](i-loopstart,j-blocks[irec_0].axis_size(1));
	
      }
    }      

  }
  
  MPI_Wait(&block_requests[isend_0],MPI_STATUS_IGNORE);
  MPI_Wait(&block_requests[jsend_1],MPI_STATUS_IGNORE);
  return;
}


void Conjugate::middle_odd_unequal(int nx)
/* only to be used when the last processor has
   less data points than the rest. */

{

  
  // write bottom left square to top right square where possible
  for (int i = 0; i < (blocks[0].axis_size(0)-1)/2; i++) {

    for (int j = 1; j < localNy/2; j++) {
      
      
      complex_update(i,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      (*ft_array)(localNz-i-blocks[1].axis_size(0)-1,localNy-j,nx) = std::conj((*ft_array)(i,j,nx));
      
    }
    /*
    nz = (globalNz- i - complex_local);
    for (int j = localNy/2 + 1; j < localNy; j++) {
      
      // since conjugatation for top, change expression for ft_noise
      ft_noise(i,j,nx) = std::conj(complex_update(nz,Ny-j,nx,ft_noise));//1.0*(Nz-i-complex_local)-1i*(1.0*(Ny-j));
      ft_noise(localNz-i-blocks[1].axis_size(0)-1,localNy-j,nx) = std::conj(ft_noise(i,j,nx));
      
      }*/
  }
  // write bottom right square to top left square where possible
  for (int i = (blocks[0].axis_size(0)-1)/2 + 1; i < blocks[0].axis_size(0); i++) {
    for (int j = 1; j < localNy/2; j++ ) {
      complex_update(i,j,nx);
      (*ft_array)(localNz-i-blocks[1].axis_size(0)-1,localNy-j,nx) = std::conj((*ft_array)(i,j,nx));
    }
  }

  
  // don't forget Nz/2 term!
  int nz = (blocks[0].axis_size(0)-1)/2;
  
  
  for (int j = 1; j < localNy/2; j++) {
    
    complex_update(nz,j,nx);//1.0*(nz+complex_local) + 1i*(1.0*j);
    (*ft_array)(nz,localNy-j,nx) = std::conj((*ft_array)(nz,j,nx));
    
  }
  
  
  // write to block and send to previous processor
  for (int i = localNz-blocks[3].axis_size(0); i < localNz; i++) {
    
    for (int j = 1; j < localNy/2; j++) {
      
      complex_update(i,j,nx);//1.0*(i+complex_local) + 1i*(1.0*j);
      blocks[3](localNz-i-1,blocks[3].axis_size(1)-j) = std::conj((*ft_array)(i,j,nx));
    }
    
  }
  
  MPI_Isend(blocks[3].data(),blocks[3].size()*2,MPI_DOUBLE,nprocs-me-2,3,
	    world,&block_requests[3]);
  
  
  // and receive extra data from the previous processor
  MPI_Recv(blocks[1].data(), blocks[1].size()*2,MPI_DOUBLE,nprocs-me-2,1,
	   world,MPI_STATUS_IGNORE);
  
  for (int i = localNz-blocks[1].axis_size(0); i < localNz; i++) {
    
    for (int j = localNy/2+1; j < localNy; j++) {
      (*ft_array)(i,j,nx) = blocks[1](i-(localNz-blocks[1].axis_size(0)),j-blocks[1].axis_size(1));
    }
  }
  
  
  MPI_Wait(&block_requests[3],MPI_STATUS_IGNORE);

  return;
  
}


void Conjugate::middle_even_unequal(int nx)
/* only to be used when the last processor has
   less data points than the rest. */

{


  // write bottom left to next processor

  for (int i = 0; i < blocks[0].axis_size(0); i++) {
    
    for (int j = 1; j < localNy/2; j++) {
      
      
      complex_update(i,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      //(ft_noise(i,j,nx) - mobility*q2*ft_nonlinear(i,j,nx)*dt + noise*dt )/(1+mobility*gamma*q2*q2);

      blocks[0](blocks[0].axis_size(0)-i-1,blocks[0].axis_size(1)-j) = std::conj((*ft_array)(i,j,nx));
      
    }

  }

  MPI_Isend(blocks[0].data(),blocks[0].size()*2,MPI_DOUBLE,nprocs-me-1,0,
	    world,&block_requests[0]);

  // write bottom left to top right in same processor where possible
  for (int i = blocks[0].axis_size(0); i < blocks[0].axis_size(0)+(blocks[1].axis_size(0)-1)/2; i++) {

    for (int j = 1; j < localNy/2; j++) {

      complex_update(i,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      (*ft_array)(localNz-(i-blocks[0].axis_size(0))-1,localNy-j,nx) = std::conj((*ft_array)(i,j,nx));

    }
  }

  // write bottom right to top left in the same processor where possible
  for (int i = blocks[0].axis_size(0)+(blocks[1].axis_size(0)-1)/2+1; i < localNz; i++) {
    for (int j = 1; j < localNy/2; j++ ) {
      complex_update(i,j,nx);
      (*ft_array)(localNz-(i-blocks[0].axis_size(0))-1,localNy-j,nx) = std::conj((*ft_array)(i,j,nx));
    }
  }
  
  // don't forget Nz/2 term (on same processor again)!
  int nz = blocks[0].axis_size(0)+(blocks[1].axis_size(0)-1)/2;
  
  
  for (int j = 1; j < localNy/2; j++) {
    
    complex_update(nz,j,nx);//1.0*(nz+complex_local) + 1i*(1.0*j);
    (*ft_array)(nz,localNy-j,nx) = std::conj((*ft_array)(nz,j,nx));
    
  }
  
  
  
  // and receive extra data from the next processor to put in top 
  MPI_Recv(blocks[2].data(), blocks[2].size()*2,MPI_DOUBLE,nprocs-me-1,2,
	   world,MPI_STATUS_IGNORE);
  
  for (int i = 0; i < blocks[2].axis_size(0); i++) {
    
    for (int j = localNy/2+1; j < localNy; j++) {
      (*ft_array)(i,j,nx) = blocks[2](i,j-blocks[2].axis_size(1));
      //      std::cout << "ft_noise(" << i+ft_array.get_local0start() << "," << j
      //		<< ") = " << (*ft_array)(i,j,nx) << std::endl;
    }
  }
  
  
  MPI_Wait(&block_requests[0],MPI_STATUS_IGNORE);

  return;
  
}




void Conjugate::line_onehalf_equal(int nx,
				   std::string half)
{

  int ny;

  int isend_0, jsend_1;
  int irec_0, jrec_1;
  
  if (half == "left") {
    isend_0 = 0;
    jsend_1 = 1;
    irec_0 = 2;
    jrec_1 = 3;
    ny = 0;
  } else if (half == "right") {
    isend_0 = 2;
    jsend_1 = 3;
    irec_0 = 0;
    jrec_1 = 1;
    ny = localNy/2;
  } else { // should never get here
    std::cout << "fucked up in coding." << std::endl;
    return;
  }

  int nz = 0;

    
  complex_update(nz,ny,nx);//1.0*(i+complex_local)+1i*(1.0*j);
  //(ft_noise(i,j,nx) - mobility*q2*ft_nonlinear(i,j,nx)*dt + noise*dt )/(1+mobility*gamma*q2*q2);
    
  lines[isend_0](0,0) = std::conj((*ft_array)(nz,ny,nx));
  
  MPI_Isend(lines[isend_0].data(),lines[isend_0].size()*2,MPI_DOUBLE,nprocs-me,isend_0,
	    world,&line_requests[isend_0]);
  
  for (int i = 1; i < localNz ; i++) {
    
    complex_update(i,ny,nx);//1.0*(i+complex_local)+1i*(1.0*j);
    lines[jsend_1](lines[jsend_1].axis_size(0)-i,0)
      = std::conj((*ft_array)(i,ny,nx));
  }
  
  MPI_Isend(lines[jsend_1].data(),lines[jsend_1].size()*2,MPI_DOUBLE,nprocs-me-1,jsend_1,
	    world,&line_requests[jsend_1]);
  
  
  
  // receive from bottom right square to fill in top left square
  // (1 <= y < Nz/2, Ny/2 < z < Ny)
  
  
  MPI_Irecv(lines[irec_0].data(), lines[irec_0].size()*2,MPI_DOUBLE,nprocs-me,irec_0,
	    world,&line_requests[irec_0]);
  
  MPI_Irecv(lines[jrec_1].data(), lines[jrec_1].size()*2,MPI_DOUBLE,nprocs-me-1,jrec_1,
	    world,&line_requests[jrec_1]);
  

  // now fill phi with lines[3] and line4 received data in an efficient manner
  
  int b2flag=0;
  int b3flag = 0;


  if (half == "left") {
    ny = localNy/2;
  } else if (half == "right") {
    ny = 0;
  } else { // should never get here
    std::cout << "fucked up in coding." << std::endl;
    return;
  }
  
  
  // next line of code ensures that whichever line data is received first
  // is filled in first
  while (!(b2flag || b3flag)) {
    MPI_Test(&line_requests[irec_0],&b2flag,MPI_STATUS_IGNORE);
    MPI_Test(&line_requests[jrec_1],&b3flag,MPI_STATUS_IGNORE);
  }
  
  
  if (b2flag) { // if lines[2] is received first

    (*ft_array)(nz,ny,nx) = lines[irec_0](nz,0);
    
    MPI_Wait(&line_requests[jrec_1],MPI_STATUS_IGNORE);
    
    for (int i = 1; i < localNz; i++) {
	(*ft_array)(i,ny,nx) = lines[jrec_1](i-1,0);
    }

  } else if (b3flag) { // if lines[3] is received first

    for (int i = 1; i < localNz; i++) {
      (*ft_array)(i,ny,nx) = lines[jrec_1](i-1,0);
    }
    
    MPI_Wait(&line_requests[irec_0],MPI_STATUS_IGNORE);

      
    (*ft_array)(nz,ny,nx) = lines[irec_0](nz,0);
  }
  
  MPI_Wait(&line_requests[isend_0],MPI_STATUSES_IGNORE); 
  MPI_Wait(&line_requests[jsend_1],MPI_STATUSES_IGNORE); 
  return;
}

void Conjugate::line_middle_odd_equal(int nx)
{

  int ny = 0;
  int nz = 0;


  // send one to next processor

  complex_update(nz,ny,nx);//1.0*(i+complex_local)+1i*(1.0*j);
  //(ft_noise(i,j,nx) - mobility*q2*ft_nonlinear(i,j,nx)*dt + noise*dt )/(1+mobility*gamma*q2*q2);
    
  lines[0](nz,0) = std::conj((*ft_array)(nz,ny,nx));

  MPI_Isend(lines[0].data(),lines[0].size()*2,MPI_DOUBLE,nprocs-me,0,
	    world,&line_requests[0]);


  // receive one from next processor


  ny = localNy/2;
  MPI_Recv(lines[2].data(), lines[2].size()*2,MPI_DOUBLE,nprocs-me,2,
	   world,MPI_STATUS_IGNORE);
    
  (*ft_array)(nz,ny,nx) = lines[2](nz,0);



  ny = 0;
  // write bottom left square to bottom right square where possible
  for (int i = 1; i < (lines[1].axis_size(0)-1)/2+1; i++) {
  
    complex_update(i,ny,nx);//1.0*(i+complex_local)+1i*(1.0*j);
    (*ft_array)(localNz-i,ny,nx) = std::conj((*ft_array)(i,ny,nx));
      
  }

  ny = localNy/2;
  // write top right square to top left square where possible
  for (int i = (lines[1].axis_size(0)-1)/2 + 2; i < localNz; i++) {

    complex_update(i,ny,nx);
    (*ft_array)(localNz-i,ny,nx) = std::conj((*ft_array)(i,ny,nx));

  }
  MPI_Wait(&line_requests[0],MPI_STATUS_IGNORE);

  return;
  
}


void Conjugate::line_middle_even_equal(int nx)
/* only to be used when all processors have same
   number of data points on them. */
{



  // don't forget Nz/2 term!

  int ny = localNy/2;

  // fill and send top right square 
    
  for (int i = 1; i < localNz; i++) {
    complex_update(i,ny,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      
    lines[3](lines[3].axis_size(0)-i,0)
	= std::conj((*ft_array)(i,ny,nx));

  }
  
  
  MPI_Isend(lines[3].data(),lines[3].size()*2,MPI_DOUBLE,nprocs-me-1,3,
	    world,&line_requests[3]);
  
  
  // receive from bottom left square, then fill in bottom right square
  
  MPI_Recv(lines[1].data(), lines[1].size()*2,MPI_DOUBLE,nprocs-me-1,1,
	   world,MPI_STATUS_IGNORE);
  
  ny = 0;
  for (int i = 1; i < localNz; i++) {
    (*ft_array)(i,ny,nx) = lines[1](i-1,0);
  }
    
  MPI_Wait(&line_requests[3],MPI_STATUS_IGNORE);
  
  return;
  
}


void Conjugate::line_onehalf_unequal( int nx,
				     std::string half)
{



  int irec_0, jrec_1, isend_0, jsend_1;
  int ny;
  
  if (half == "left") {
    isend_0 = 0;
    jsend_1 = 1;
    irec_0 = 2;
    jrec_1 = 3;
    ny = 0;
  } else if (half == "right") {
    isend_0 = 2;
    jsend_1 = 3;
    irec_0 = 0;
    jrec_1 = 1;
    ny = localNy/2;
  } else {
    std::cout << "fucked up somehwere" << std::endl;
    return;
  }


  int loopstart = 0;
  if (me ==0 ) {
    loopstart = 1;
  }

  for (int i = loopstart; i < lines[isend_0].axis_size(0)+loopstart; i++) {

    complex_update(i,ny,nx);//1.0*(i+complex_local)+1i*(1.0*j);
    //(ft_noise(i,j,nx) - mobility*q2*ft_nonlinear(i,j,nx)*dt + noise*dt )/(1+mobility*gamma*q2*q2);

    lines[isend_0](lines[isend_0].axis_size(0)-i-(1-loopstart),0) = std::conj((*ft_array)(i,ny,nx));
      
  }
  
  MPI_Isend(lines[isend_0].data(),lines[isend_0].size()*2,MPI_DOUBLE,nprocs-me-1,isend_0,
	    world,&line_requests[isend_0]);
  
  for (int i = lines[isend_0].axis_size(0)+loopstart; i < localNz ; i++) {

    complex_update(i,ny,nx);//1.0*(i+complex_local)+1i*(1.0*j);
    lines[jsend_1](lines[jsend_1].axis_size(0)-(i-lines[isend_0].axis_size(0)+1-loopstart),0)
      = std::conj((*ft_array)(i,ny,nx));
  }
  
  MPI_Isend(lines[jsend_1].data(),lines[jsend_1].size()*2,MPI_DOUBLE,nprocs-me-2,jsend_1,
	    world,&line_requests[jsend_1]);
  
  
  
  // receive from top right side line to fill in top left line
  // (1 <= y < Nz/2, z = Ny/2)
  

  if (half == "left") {
    ny = localNy/2;
  } else if (half == "right") {
    ny = 0;
  } else { // should never get here
    std::cout << "fucked up in coding." << std::endl;
    return;
  }

  
  MPI_Irecv(lines[irec_0].data(), lines[irec_0].size()*2,MPI_DOUBLE,nprocs-me-1,irec_0,
	    world,&line_requests[irec_0]);
  
  MPI_Irecv(lines[jrec_1].data(), lines[jrec_1].size()*2,MPI_DOUBLE,nprocs-me-2,jrec_1,
	    world,&line_requests[jrec_1]);
  
  
  // now fill phi with lines[3] and line4 received data in an efficient manner
  
  int b2flag=0;
  int b3flag = 0;
  
  
  // next line of code ensures that whichever line data is received first
  // is filled in first
  while (!(b2flag || b3flag)) {
    MPI_Test(&line_requests[irec_0],&b2flag,MPI_STATUS_IGNORE);
    MPI_Test(&line_requests[jrec_1],&b3flag,MPI_STATUS_IGNORE);
  }
  
  
  if (b2flag) { // if lines[irec_0] is received first
    for (int i = loopstart; i < lines[irec_0].axis_size(0)+loopstart; i++) {
      (*ft_array)(i,ny,nx) = lines[irec_0](i-loopstart,0);
      
    }

    
    MPI_Wait(&line_requests[jrec_1],MPI_STATUS_IGNORE);
    
    for (int i = lines[irec_0].axis_size(0)+loopstart; i < localNz; i++) {
      (*ft_array)(i,ny,nx) = lines[jrec_1](i-(lines[irec_0].axis_size(0)+loopstart),0);
    }
  } else if (b3flag) { // if lines[3] is received first

    for (int i = lines[irec_0].axis_size(0)+loopstart; i < localNz; i++) {
	(*ft_array)(i,ny,nx) = lines[jrec_1](i-(lines[irec_0].axis_size(0)+loopstart),0);
    }
    
    MPI_Wait(&line_requests[irec_0],MPI_STATUS_IGNORE);
    
    
    for (int i = loopstart; i < lines[irec_0].axis_size(0)+loopstart; i++) {
      (*ft_array)(i,ny,nx) = lines[irec_0](i-loopstart,0);
	
    }

  }
  
  MPI_Wait(&line_requests[irec_0],MPI_STATUS_IGNORE);    
  MPI_Wait(&line_requests[jrec_1],MPI_STATUS_IGNORE);

  return;
  

}

void Conjugate::line_middle_odd_unequal(int nx)
{

  int ny = 0;
  
  // write bottom left square to bottom right square where possible
  for (int i = 0; i < (lines[0].axis_size(0)-1)/2; i++) {
    
    complex_update(i,ny,nx);//1.0*(i+complex_local)+1i*(1.0*j);
    (*ft_array)(localNz-i-lines[1].axis_size(0)-1,ny,nx) = std::conj((*ft_array)(i,ny,nx));
      
  }

  // write top right square to top left square where possible

  ny = localNy/2;
  for (int i = (lines[0].axis_size(0)-1)/2 + 1; i < lines[0].axis_size(0); i++) {
    complex_update(i,ny,nx);
    (*ft_array)(localNz-i-lines[1].axis_size(0)-1,ny,nx) = std::conj((*ft_array)(i,ny,nx));
  }

  
  // write to line and send to previous processor

  ny = localNy/2;
  
  for (int i = localNz-lines[3].axis_size(0); i < localNz; i++) {
    complex_update(i,ny,nx);//1.0*(i+complex_local) + 1i*(1.0*j);
    lines[3](localNz-i-1,0) = std::conj((*ft_array)(i,ny,nx));
  }
  
  MPI_Isend(lines[3].data(),lines[3].size()*2,MPI_DOUBLE,nprocs-me-2,3,
	    world,&line_requests[3]);
  
  
  // and receive extra data from the previous processor
  MPI_Recv(lines[1].data(), lines[1].size()*2,MPI_DOUBLE,nprocs-me-2,1,
	   world,MPI_STATUS_IGNORE);

  ny = 0;
  for (int i = localNz-lines[1].axis_size(0); i < localNz; i++) {
    (*ft_array)(i,ny,nx) = lines[1](i-(localNz-lines[1].axis_size(0)),0);
  }
  
  
  MPI_Wait(&line_requests[3],MPI_STATUS_IGNORE);

  return;
  
}

void Conjugate::line_middle_even_unequal(int nx)
{


  int ny = 0;

  // write bottom left and send to next processor
  for (int i = 0; i < lines[0].axis_size(0); i++) {
      
    complex_update(i,ny,nx);//1.0*(i+complex_local)+1i*(1.0*j);


    lines[0](lines[0].axis_size(0)-i-1,0) = std::conj((*ft_array)(i,ny,nx));
      
  }

  // 

  MPI_Isend(lines[0].data(),lines[0].size()*2,MPI_DOUBLE,nprocs-me-1,0,
	    world,&line_requests[0]);

  // write bottom left to bottom right in same processor where possible
  for (int i = lines[0].axis_size(0); i < lines[0].axis_size(0)+(lines[1].axis_size(0)-1)/2; i++) {
    complex_update(i,ny,nx);//1.0*(i+complex_local)+1i*(1.0*j);
    (*ft_array)(localNz-(i-lines[0].axis_size(0))-1,ny,nx) = std::conj((*ft_array)(i,ny,nx));
  }

  ny = localNy/2;
  // write top right to top left in the same processor where possible
  for (int i = lines[0].axis_size(0)+(lines[1].axis_size(0)-1)/2+1; i < localNz; i++) {
    complex_update(i,ny,nx);
    (*ft_array)(localNz-(i-lines[0].axis_size(0))-1,ny,nx) = std::conj((*ft_array)(i,ny,nx));
  }

  // and receive extra data from the next processor to put in left

  ny = localNy/2;
  MPI_Recv(lines[2].data(), lines[2].size()*2,MPI_DOUBLE,nprocs-me-1,2,
	   world,MPI_STATUS_IGNORE);
  
  for (int i = 0; i < lines[2].axis_size(0); i++) {
    (*ft_array)(i,ny,nx) = lines[2](i,0);
  }

  MPI_Wait(&line_requests[0],MPI_STATUS_IGNORE);


  return;
  
}


