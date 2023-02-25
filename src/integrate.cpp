#include "integrate.hpp"

#include "domain.hpp"
#include "grid.hpp"
#include "comm_brick.hpp"
#include "atom.hpp"
#include "neighbor.hpp"
#include "fix.hpp"
#include "pair.hpp"
#include "compute.hpp"
#include "dump.hpp"
#include "ps_pde/fftw_mpi_3darray.hpp"

#include <mpi.h>
#include <iostream>
#include <string>
#include <set>
#include <chrono>

using namespace PHAFD_NS;

Integrate::Integrate(PHAFD *phafd) : Pointers(phafd) {
  nsteps = 0;
  firststep = 0;
  timestep = 0;
  dt = 0.0;
  
}

Integrate::~Integrate() = default;

void Integrate::setup()
{  

  // need to apply pbc to atoms in case they are outside of simulation domain
  domain->pbc();
  // can only call borders after atoms are in the simulation domain
  commbrick->borders();

  // can only call neighbors after pair coefficients have been set (for force cutoff)
  neighbor->build();

  
  for (auto &dump : dumps) {
    dump->setup();
    dump->write_collection_header();
  }

  
  for (auto &fix : fixes)
    fix->reset_dt();

  for (auto &fix : fixes)
    fix->setup();


  // do the below to potentially save things on step zero.

  for (auto &compute: computes)
    compute->start_of_step();
  
  for (auto &fix: fixes)
    fix->start_of_step();
  
  
  for (auto &dump : dumps)
    dump->start_of_step();
  
  
  
  atoms->Fs.setZero();
  grid->nonlinear->setZero();
  
  for (auto & pair : pairs)
    pair->compute();
  
  commbrick->reverse_comm();
  
  // additional terms from e.g. chemical potential, which are added to grid->nonlinear
  for (auto &fix: fixes)
    fix->post_force();


  // mainly just fourier transforming phi and nonlinear
  for (auto &fix : fixes)
    fix->pre_final_integrate();


  // computes which act on fourier space grids
  for (auto &compute : computes)
    compute->in_fourier();


  // mainly just inverse fourier transforming
  for (auto &fix : fixes)
    fix->post_final_integrate();

  for (auto &compute : computes)
    compute->end_of_step();

  for (auto &fix : fixes)
    fix->end_of_step();
    
    
  for (auto &dump :dumps)
    if (timestep % dump->every == 0) {
      if (commbrick->me == 0)
	std::cout << "Saving on step " << timestep << std::endl;
      
      
      dump->write_collection_middle();
      
      
    }
  

}

void Integrate::run()
{

  if (commbrick->me == 0) 
    std::cout << "Running simulation of solution." << std::endl;
  
    
  
  int errflag = 0;
  int total_errflag = 0;

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();



  for (int i = 0; i < nsteps; i++) {
  
    timestep ++;

    
    for (auto &compute: computes)
      compute->start_of_step();
    
    for (auto &fix: fixes)
      fix->start_of_step();


    for (auto &dump : dumps)
      dump->start_of_step();

    
    if (neighbor->decide()) {
      domain->pbc();
      commbrick->borders();
      neighbor->build();
    } else {
      commbrick->forward_comm();
    }


    atoms->Fs.setZero();

    for (auto & pair : pairs)
      pair->compute();

    commbrick->reverse_comm();

    for (auto &fix : fixes)
      fix->initial_integrate();

    grid->nonlinear->setZero();
    atoms->Fs.setZero();
    
    for (auto & pair : pairs)
      pair->compute();

    commbrick->reverse_comm();

    // additional terms from e.g. chemical potential, which are added to grid->nonlinear
    for (auto &fix: fixes)
      fix->post_force();


    // mainly just fourier transforming phi and nonlinear
    for (auto &fix : fixes)
      fix->pre_final_integrate();


    // computes which act on fourier space grids
    for (auto &compute : computes)
      compute->in_fourier();

    
    // updating stuff (including fourier transformed phi and nonlinear
    for (auto &fix : fixes)
      fix->final_integrate();


    // mainly just inverse fourier transforming
    for (auto &fix : fixes)
      fix->post_final_integrate();

    for (auto &compute : computes)
      compute->end_of_step();
    
    for (auto &fix : fixes)
      fix->end_of_step();
    
    
    for (auto &dump :dumps)
      if (timestep % dump->every == 0) {
	if (commbrick->me == 0)
	  std::cout << "Saving on step " << timestep << std::endl;


	dump->write_collection_middle();

      
      }
    
    if (std::isnan((*grid->phi)(0,0,0)))
      errflag = 1;
	
    MPI_Allreduce(&errflag, &total_errflag,1,MPI_INT,MPI_SUM,world);

    if (total_errflag) {

      for (auto &dump : dumps)
	dump->write_collection_footer();
      
      throw std::runtime_error("NAN encountered in phi.");
    }    
    
  }

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

  for (auto &dump :dumps)
    dump->write_collection_footer();

  finalise();
  
  if (commbrick->me == 0) {
    std::cout << "Run time = "
	      << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1e6
	      << "seconds." << std::endl;  
  }


}

void Integrate::run_until_touching(double touch_point)
{

  if (commbrick->me == 0) 
    std::cout << "Running simulation of solution." << std::endl;
  

  for (auto & compute : computes)
    if (compute->per_ftgrid)
      throw std::runtime_error("cannot compute fourier properties for run_until_touching.");

  for (auto & fix : fixes)
    if (fix->averaging)
      throw std::runtime_error("cannot use fix/ave for run_until_touching.");
  
  int errflag = 0;
  int total_errflag = 0;

  int breakflag = 0;
  int total_breakflag = 0;

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();



  for (int i = 0; i < nsteps; i++) {
  
    timestep ++;

    
    for (auto &compute: computes)
      compute->start_of_step();
    
    for (auto &fix: fixes)
      fix->start_of_step();


    for (auto &dump : dumps)
      dump->start_of_step();

    
    if (neighbor->decide()) {
      domain->pbc();
      commbrick->borders();
      neighbor->build();
    } else {
      commbrick->forward_comm();
    }


    atoms->Fs.setZero();

    for (auto & pair : pairs)
      pair->compute();

    commbrick->reverse_comm();

    for (auto &fix : fixes)
      fix->initial_integrate();

    grid->nonlinear->setZero();
    atoms->Fs.setZero();
    
    for (auto & pair : pairs)
      pair->compute();

    commbrick->reverse_comm();

    // additional terms from e.g. chemical potential, which are added to grid->nonlinear
    for (auto &fix: fixes)
      fix->post_force();


    // mainly just fourier transforming phi and nonlinear
    for (auto &fix : fixes)
      fix->pre_final_integrate();


    // computes which act on fourier space grids
    for (auto &compute : computes)
      compute->in_fourier();

    
    // updating stuff (including fourier transformed phi and nonlinear
    for (auto &fix : fixes)
      fix->final_integrate();


    // mainly just inverse fourier transforming
    for (auto &fix : fixes)
      fix->post_final_integrate();

    for (auto &compute : computes)
      compute->end_of_step();
    
    for (auto &fix : fixes)
      fix->end_of_step();
    
    
    for (auto &dump :dumps)
      if (timestep % dump->every == 0) {
	if (commbrick->me == 0)
	  std::cout << "Saving on step " << timestep << std::endl;


	dump->write_collection_middle();

      
      }
    
    if (std::isnan((*grid->phi)(0,0,0)))
      errflag = 1;
	
    MPI_Allreduce(&errflag, &total_errflag,1,MPI_INT,MPI_SUM,world);

    if (total_errflag) {

      for (auto &dump : dumps)
	dump->write_collection_footer();
      
      throw std::runtime_error("NAN encountered in phi.");
    }


    
    if (commbrick->me == 0) {

      if ((atoms->xs.col(atoms->nowned-1)-atoms->xs.col(0)).norm() < touch_point) {
	std::cout << "cyclisation occured at step "
		  << timestep << std::endl;

	breakflag = 1;
      }

    }

    MPI_Allreduce(&breakflag, &total_breakflag,1,MPI_INT,MPI_SUM,world);

    if (total_breakflag) {

      for (auto &dump : dumps) 
	dump->require_calculations();

      for (auto &compute : computes)
	compute->end_of_step();
    
      for (auto &fix : fixes)
	fix->end_of_step();

      for (auto &dump : dumps)
	dump->write_collection_middle();
	

      break;

    }
    
    
    
  }

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

  for (auto &dump :dumps)
    dump->write_collection_footer();
  
  if (commbrick->me == 0) {
    std::cout << "number of neighbor calls = " << neighbor->get_ncalls() << std::endl;
    std::cout << "other flags:" << std::endl;
  }


  finalise();
  
  if (commbrick->me == 0) {  
    std::cout << "Run time = "
	      << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1e6
	      << "seconds." << std::endl;  
  }


}


void Integrate::finalise()
{
  std::vector<int> num_less_total,num_great_total;

  for (auto &fix : fixes)  {
    int tmp;
    MPI_Reduce(&(fix->num_less_zero),&tmp,1,MPI_INT,MPI_SUM,0,world);
    num_less_total.push_back(tmp);
    MPI_Reduce(&(fix->num_great_zero),&tmp,1,MPI_INT,MPI_SUM,0,world);
    num_great_total.push_back(tmp);
  }

  
  if (commbrick->me == 0) {
    
    int count = 0;
    for (auto &fix : fixes)  {
      std::cout << fix->name << " had phi < 0 " << num_less_total.at(count) << " times, and "
		<< "phi > 1 " << num_great_total.at(count) << "times." << std::endl;
      count += 1;
    }
    
  }

  return;
}
