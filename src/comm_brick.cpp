#include <mpi.h>

#include "atom.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "neighbor.hpp"

#define BIG 1e20


/*
  Assumes that the processors cover all of x y space, but are partitioned along
  the z axis at regular intervals.
  
  rebuild neighbor lists from the assumption that each processor owns npoly_id
  polymers, which are stored as the first nowned atoms in the atom class. This
  means that a processor has two distinct types of atoms: those which are part
  of one of the npoly_id polymers (total nowned atoms) and those which are
  not part of the npoly_id atoms, but are in the local spatial region of the
  processor.
  
  So, an owned atom of a processor is shared with any other processor for which
  the atom is spatially close enough (i.e. within a specified cutoff distance).

  So, each processor has nprocs-1 swaps, one for each of the other processors.
  It iterates through these swaps to receive and send atoms to the other
  processors when correct to do so.
  
*/


using namespace PHAFD_NS;

CommBrick::CommBrick(PHAFD *phafd)
  : Pointers(phafd), xyswaps(4)
{

  MPI_Comm_size(world,&nprocs);
  MPI_Comm_rank(world,&me);
  
  maxneed[0] = 2;
  maxneed[1] = 2;
  maxneed[2] = nprocs;
  
  // set values as if there were no processors at all

  zprd = 0.0;
  zinterval = 0.0;
  zlo = 0.0;
  cutghost = 0.0;

  size_forward=0;
  size_reverse=0;
  size_border=0;
  


  // there should be 2 swaps in the x dim, 2 swaps in the y dim, and
  // nprocs swaps in the z dim

  
  pbc.resize(xyswaps);
  
  nswap = nprocs + xyswaps;
  slablo.resize(nswap);
  slabhi.resize(nswap);


    

  subloz.resize(nprocs);
  subhiz.resize(nprocs);

  
  sendproc.resize(nswap);
  recvproc.resize(nswap);
  firstrecv.resize(nswap);
  sendnum.resize(nswap);
  recvnum.resize(nswap);

  size_reverse_recv.resize(nswap);
  size_forward_recv.resize(nswap);
  size_reverse_send.resize(nswap);

  sendlists.resize(nswap);

  pbcz.resize(nprocs);  

  
  
}

CommBrick::~CommBrick() = default;

/* must be called after domain subbox is created, and after neighbor setup is called. */
void CommBrick::setup()
{

  const int zdim = 2;

  cutghost = neighbor->cutneigh;
  zprd = domain->period[2];
  zinterval = domain->subhi[2]-domain->sublo[2];
  zlo = domain->boxlo[2];


  // start with zdim, which uses different structure than the x and y dim swaps
  // (the x and y dim swaps use the same method as LAMMPS)

  int iswap = 0;

  // z dimension




  for (int ineed = 0; ineed < maxneed[zdim]; ineed++ ) {
 
    sendproc[iswap] = (me - 1 - iswap < 0 ? me + nprocs -  iswap - 1  : me - 1  - iswap );
    recvproc[iswap] = (me + 1 + iswap) % nprocs;
    subloz[iswap] = sendproc[iswap]*zinterval + zlo;
    subhiz[iswap] = (sendproc[iswap]+1)*zinterval + zlo;

    
    slablo[iswap] = (sendproc[iswap] == 0 ? nprocs : sendproc[iswap] )*zinterval
      - cutghost + zlo;
    slabhi[iswap] = ((sendproc[iswap]+1) % nprocs)*zinterval + cutghost
      + zlo;
     
    
    iswap ++;
    
  }


  
  // x and y dimensions, only one processor in each
  
  for (int dim = 1; dim >= 0; dim--) {

    pbc[iswap-nprocs][0] = pbc[iswap-nprocs][1] = pbc[iswap-nprocs][2] = 0.0;
    
    for (int ineed = 0; ineed < maxneed[dim]; ineed++ ) {

	
      sendproc[iswap] = me;
      recvproc[iswap] = me;
      

      if (ineed % 2 == 0) {
	slablo[iswap] =  -BIG;
	slabhi[iswap] = domain->boxlo[dim] + cutghost;
	pbc[iswap-nprocs][dim] = domain->period[dim];
      }
      else {
	slablo[iswap] = domain->boxhi[dim] - cutghost;
	slabhi[iswap] = BIG;
	pbc[iswap-nprocs][dim] = -domain->period[dim];
      }


      
      iswap ++;
    }
  }

  return;
}


void CommBrick::forward_comm()
{

  double *buf;
  MPI_Request request;

  int n;
  
  for (int iswap = 0; iswap < nswap; iswap++) {
    if (sendproc[iswap] != me) {
      if (size_forward_recv[iswap]) {
	buf = &(atoms->xs(0,firstrecv[iswap]));
	MPI_Irecv(buf,size_forward_recv[iswap],MPI_DOUBLE,recvproc[iswap],0,world,&request);
      }
      buf_send.resize(sendnum[iswap]*size_forward);
      if (iswap < nprocs)
	n = atoms->pack_comm_in_z(sendnum[iswap],sendlists[iswap],buf_send.data(),pbcz[iswap]);
      else
	n = atoms->pack_comm(sendnum[iswap],sendlists[iswap],buf_send.data(),pbc[iswap-nprocs]);
      if (n) MPI_Send(buf_send.data(),n,MPI_DOUBLE,sendproc[iswap],0,world);
      if (size_forward_recv[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
    } else {
      if (sendnum[iswap]) {
	if (iswap < nprocs)
	  n = atoms->pack_comm_in_z(sendnum[iswap],sendlists[iswap],
				   &(atoms->xs(0,firstrecv[iswap])),pbcz[iswap]);
	else
	  n = atoms->pack_comm(sendnum[iswap],sendlists[iswap],
			      &(atoms->xs(0,firstrecv[iswap])),pbc[iswap-nprocs]);
      }
    }
  }

  return;
}



void CommBrick::reverse_comm()
{
  int n;
  MPI_Request request;
  double *buf;

  // exchange data with another proc
  // if other proc is self, just copy

  for (int iswap = nswap-1; iswap >= 0; iswap--) {
    if (sendproc[iswap] != me) {
      
      if (size_reverse_recv[iswap]) {
	buf_recv.resize(size_reverse_recv[iswap]*size_reverse);
	MPI_Irecv(buf_recv.data(),size_reverse_recv[iswap],MPI_DOUBLE,sendproc[iswap],
		  0,world,&request);
      }
      if (size_reverse_send[iswap]) {
	buf = &(atoms->Fs(0,firstrecv[iswap]));
	MPI_Send(buf,size_reverse_send[iswap],MPI_DOUBLE,recvproc[iswap],0,world);
      }
      if (size_reverse_recv[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
      atoms->unpack_reverse(sendnum[iswap],sendlists[iswap],buf_recv.data());
      
    } else {
      if (sendnum[iswap])
	atoms->unpack_reverse(sendnum[iswap],sendlists[iswap],&(atoms->Fs(0,firstrecv[iswap])));
      
    }
  }
}




void CommBrick::borders()
{

  size_forward = atoms->get_size_fwd();
  size_reverse = atoms->get_size_rev();
  size_border = atoms->get_size_bdr();
  
  int iswap = 0;

  MPI_Request request;

  int nsend,nrecv;
  int nlocalsend,nlocalrecv;
  
  int nfirst,nlast;


  atoms->nlocal = atoms->ngathered = atoms->nghost = 0;


  for (int dim = 2; dim >= 0; dim -- ) {


    nlast = atoms->nowned;

    for (int ineed = 0; ineed < maxneed[dim]; ineed++) {

      double lo = slablo[iswap];
      double hi = slabhi[iswap];

      sendlists[iswap].clear();
      nlocalsend = nlocalrecv = 0;

      if (dim == 2)  {
	labels.clear();      
	pbcz[iswap].clear();
	if (sendproc[iswap] == 0) {
	  // need this for periodic boundary conditions
	  for (int i = 0; i < atoms->nowned; i++) {
	    if (atoms->xs(dim,i) <= hi) {
	      sendlists[iswap].push_back(i);
	      pbcz[iswap].push_back(0);
	      labels.push_back(Atom::GHOST);
	      if (atoms->xs(dim,i) <= subhiz[iswap]) {
		nlocalsend ++;
		labels.back() = Atom::LOCAL;
	      }
	    } else if (atoms->xs(dim,i) >= lo) {
	      sendlists[iswap].push_back(i);
	      pbcz[iswap].push_back(-zprd);
	      labels.push_back(Atom::GHOST);
	    }
	  }
	}  else if (sendproc[iswap] == nprocs - 1) {
	  // need this for periodic boundary conditions
	  for (int i = 0; i < atoms->nowned; i++) {
	    if (atoms->xs(dim,i) >= lo) {
	      sendlists[iswap].push_back(i);
	      pbcz[iswap].push_back(0);
	      labels.push_back(Atom::GHOST);
	      if (atoms->xs(dim,i) >= subloz[iswap]) {
		nlocalsend ++;
		labels.back() = Atom::LOCAL;
	      }
	    } else if (atoms->xs(dim,i) <= hi) {
	      sendlists[iswap].push_back(i);
	      pbcz[iswap].push_back(zprd);
	      labels.push_back(Atom::GHOST);
	    }
	  }
	} else {
	  for (int i = 0; i < atoms->nowned; i++) {
	    if (atoms->xs(dim,i) >= lo && atoms->xs(dim,i) <= hi) {
	      sendlists[iswap].push_back(i);
	      pbcz[iswap].push_back(0);
	      labels.push_back(Atom::GHOST);
	      if (atoms->xs(dim,i) >= subloz[iswap] &&
		  atoms->xs(dim,i) <= subhiz[iswap]) {
		nlocalsend ++;
		labels.back() = Atom::LOCAL;
	      }
	    }
	  }
	}
      } else {
	if (ineed % 2 == 0) {
	  nfirst = nlast;
	  nlast = atoms->nowned + atoms->ngathered;
	}

	for (int i = nfirst; i < nlast; i++) {
	  if (atoms->xs(dim,i) >= lo && atoms->xs(dim,i) <= hi) 
	    sendlists[iswap].push_back(i);
	}

      }
      nsend = sendlists[iswap].size();
      buf_send.resize(nsend*size_border);


      
      int n;

      if (dim == 2) 
	n = atoms->pack_border_in_z(nsend,sendlists[iswap],buf_send.data(),pbcz[iswap],
				   labels);
      else
	n = atoms->pack_border(nsend,sendlists[iswap],buf_send.data(),pbc[iswap-nprocs]);

      if (sendproc[iswap] != me) {
	MPI_Sendrecv(&nsend,1,MPI_INT,sendproc[iswap],0,
		     &nrecv,1,MPI_INT,recvproc[iswap],0,world,MPI_STATUS_IGNORE);
	if (dim == 2) {
	  MPI_Sendrecv(&nlocalsend,1,MPI_INT,sendproc[iswap],0,
		       &nlocalrecv,1,MPI_INT,recvproc[iswap],0,world,MPI_STATUS_IGNORE);
	}

	buf_recv.resize(nrecv*size_border);
	if (nrecv) MPI_Irecv(buf_recv.data(),buf_recv.size(),MPI_DOUBLE,
			     recvproc[iswap],0,world,&request);
	if (n) MPI_Send(buf_send.data(),n,MPI_DOUBLE,sendproc[iswap],0,world);
	if (nrecv) MPI_Wait(&request,MPI_STATUS_IGNORE);
	atoms->unpack_border(nrecv,atoms->nowned+atoms->ngathered,buf_recv.data());
      } else {
	nrecv = nsend;
	nlocalrecv = nlocalsend;
	atoms->unpack_border(nrecv,atoms->nowned+atoms->ngathered,buf_send.data());
      }

      sendnum[iswap] = nsend;
      recvnum[iswap] = nrecv;
      size_forward_recv[iswap] = nrecv*size_forward;
      size_reverse_send[iswap] = nrecv*size_reverse;
      size_reverse_recv[iswap] = nsend*size_reverse;
      firstrecv[iswap] = atoms->nowned+atoms->ngathered;
      atoms->ngathered += nrecv;
      atoms->nlocal += nlocalrecv;

      iswap ++;
    }	  
  }

  atoms->nghost = atoms->ngathered-atoms->nlocal;
  return;
      
}
