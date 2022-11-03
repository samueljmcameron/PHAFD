#include <iostream>
#include <map>
#include <complex>
#include <chrono>
#include <cstring>

#include <fstream>

#include <fftw3-mpi.h>
#include "ps_pde/solutionparams.hpp"
#include "ps_pde/input.hpp"
#include "beadrodpmer/input.hpp"

#include "run.hpp"

int main(int argc, char **argv)
{

  

  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Info info = MPI_INFO_NULL;

  int ierr = MPI_Init(NULL,NULL);
  
  if (ierr != MPI_SUCCESS) {
    std::cout << "Error starting MPI program. Terminating." << std::endl;
    MPI_Abort(MPI_COMM_WORLD,ierr);
  }
  int mpi_size, id;

  ierr = MPI_Comm_size(comm,&mpi_size);
  ierr = MPI_Comm_rank(comm,&id);

  fftw_mpi_init();

  std::ifstream infile,nucfile;

  std::map<std::string,std::string> variables;

  std::string simulation_type = "run";

  int iarg = 1;  
  while(iarg < argc) {
    if (strcmp(argv[iarg],"-in") == 0) {
      if (iarg+1 == argc) {
	std::cerr << "Error: input flag '-in' specified, but no file given."
		  << std::endl;
	return EXIT_FAILURE;
      }
      infile.open(argv[iarg+1]);
      iarg += 2;
      
    } else if (strcmp(argv[iarg],"-nuc") == 0) {
      if (iarg+1 == argc) {
	std::cerr << "Error: input flag '-nuc' specified, but no file given."
		  << std::endl;
	return EXIT_FAILURE;
      }
      nucfile.open(argv[iarg+1]);
      iarg += 2;
      
    } else if (strcmp(argv[iarg],"-var") == 0) {
      
      if (iarg + 2 >= argc) {
	std::cerr << "Error: invalid command line variable specification."
		  << std::endl;
	return EXIT_FAILURE;
      }
      variables[argv[iarg+1]] = argv[iarg+2];
      iarg += 3;
    } else {
      std::cerr << "Error: invalid command line variable specification."
		<< std::endl;
      return EXIT_FAILURE;
    }
  }

  
  if (not infile.is_open()) {
    std::cerr << "Error: need to specify input file." << std::endl;
    return EXIT_FAILURE;
  }

  

  std::string line;
  std::vector<std::string> splitvec;


  std::vector<std::vector<std::string>> polymersplitvecs;

    
  // store global parameters
  psPDE::GlobalParams gp(comm,id,mpi_size,infile,variables,line,
			 {"build_solution","no_tether","single_tether","double_tether",
			     "read_no_tether","read_single_tether",
			     "read_double_tether"});
  
  // but still need to get the local parameters for both the 


  std::vector<std::string> polymertypes;

  psPDE::SolutionParams solparams;

  do {

    line = line.substr(0,line.find_first_of("#"));
    splitvec = psPDE::input::split_line(line);
    

    

    if (splitvec.size() == 0) {
      continue;
    } else if (splitvec[0] == "build_solution") {


      splitvec.erase(splitvec.begin());
      for (auto &c : splitvec) { 
	psPDE::input::convertVariable(c,variables);

      }


    
      solparams = psPDE::SolutionParams(splitvec);
      
      
    } else if (splitvec[0] == "single_tether" ||
	       splitvec[0] == "double_tether" ||
	       splitvec[0] == "no_tether") {

      
      polymertypes.push_back(splitvec.at(0));
      
	  
      splitvec.erase(splitvec.begin());
      
      for (auto &c : splitvec) {
      	BeadRodPmer::input::convertVariable(c,variables);
	
      }

      polymersplitvecs.push_back(splitvec);


    } else if (splitvec[0] == "read_single_tether" ||
	       splitvec[0] == "read_double_tether" ||
	       splitvec[0] == "read_no_tether") {

      if (!gp.read_flag) 
	throw std::runtime_error("Cannot read polymer when not reading in data from" 
				 " previous simulation.");



      
      for (auto &c : splitvec) {
      	BeadRodPmer::input::convertVariable(c,variables);
      }

      polymertypes.push_back(splitvec.at(0)+splitvec.at(1));
      
      splitvec.erase(splitvec.begin());
      splitvec.erase(splitvec.begin());


      polymersplitvecs.push_back(splitvec);

      
    }
      
      
  } while (std::getline(infile,line)); 

  if (id == 0) {
    gp.printall();
    solparams.printall();
  }

  std::vector<std::vector<double>> X_is;
  

  std::vector<double> radii;
  std::vector<double> viscosities;
  
  if (not nucfile.is_open() && not gp.read_flag && not gp.restart_flag) {
    std::cerr << "Warning! No nucleation file specified, so no nucleation sites will be generated."
	      << std::endl;
  } else {

    double X_x,X_y,X_z;
    double viscosity,radius;
    while(nucfile >> X_x >> X_y >> X_z >> radius >> viscosity) {
      X_is.push_back({X_x,X_y,X_z});
      radii.push_back(radius);
      viscosities.push_back(viscosity);      
    }

  }

  // now, give the solution parameters the list of polymers to include, and the
  // nucleation sites, run the simulation

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  run(gp,solparams,polymertypes,polymersplitvecs,X_is,radii,viscosities);

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Run time = "
	    << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1e6
	    << "seconds." << std::endl;  
  
  fftw_mpi_cleanup();

  ierr = MPI_Finalize();
  
  return 0;
}
