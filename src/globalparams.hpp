#ifndef GLOBALPARAMS_HPP
#define GLOBALPARAMS_HPP

#include <fstream>
#include <set>
#include <map>
#include <mpi.h>


#include "ps_pde/griddata.hpp"

class GlobalParams {
private:

  std::set<std::string> pset {"steps","volFrac", "dt", "solution_dump_every",
      "solution_dump_file", "thermo_every", "polymer_dump_every",
      "polymer_dump_file", "thermo_file","boxgrid", "boxdims",
      "solution_seed","restart","polymer_equilibration"};

  const int default_steps;
  const int default_grid;
  const double default_length;
  const double default_dt;
  const int default_dump_every;
  const int default_thermo_every;
  const std::string default_solution_dump_file;
  const std::string default_polymer_dump_file;
  const std::string default_thermo_file;
  const double default_volFrac;
  const double default_variance;
  const int default_seed;
  const int default_polymer_equilibration;


  bool fft_is_transposed;
public:
  GlobalParams(const MPI_Comm, const int, const int,std::ifstream&,
	       std::map<std::string,std::string> const &,
	       std::string&);

  bool restart_flag;
  int startstep;
  double starttime;

  
  int steps,seed,solution_dump_every,polymer_dump_every,thermo_every;
  int polymer_equilibration;
  double dt, volFrac,variance;
  std::string solution_dump_file,polymer_dump_file,thermo_file;
  
  const MPI_Comm comm;
  const int id, mpi_size;

  psPDE::GridData realspace;
  psPDE::GridData fourier;
  


  void printall();

  
  
};


#endif
