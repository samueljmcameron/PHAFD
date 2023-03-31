#include <iostream>
#include <algorithm>
#include <set>

#include "utility.hpp"
#include "atom.hpp"

#include "read_atoms.hpp"
 



using namespace PHAFD_NS;

ReadAtoms::ReadAtoms(PHAFD *phafd) : Pointers(phafd) {};

int ReadAtoms::read_file(const std::string & fname_in)
{


  fname = fname_in;
  int errflag = SUCCESS;
  int total_errflag;

  // open the file to be read
  
  datafile.open(fname);
  
  if (datafile.fail()) {
    std::cerr << "could not open file " << fname << std::endl;
    errflag = NOFILE;
  }

  if (errflag != SUCCESS) errflag = 1;
  else errflag = 0;

  MPI_Allreduce(&errflag,&total_errflag,1,MPI_INT,MPI_MAX,world);
  
  if (total_errflag) return NOFILE;

  // create the atoms from the file info
  
  errflag = create_atoms(atoms->xs.cols());
  
  MPI_Allreduce(&errflag,&total_errflag,1,MPI_INT,MPI_SUM,world);

  if (total_errflag) return FORMAT_ERROR;
  
  datafile.close();
  datafile.clear();
  
  return SUCCESS;
  
}



int ReadAtoms::create_atoms(int totalatoms)
/* Construct atoms from the lines of the data file. */
{

  // allocate unique_ptr for atoms



  int iatom = 0;


  std::string line;

  std::vector<std::string> v_line;

  int created_atoms;


  while (std::getline(datafile,line)) {
      
    if (line == "" || line[0] == '#') continue;

    v_line = utility::split_line(line);
    
	  
    created_atoms = atoms->add_atom(v_line,iatom);
	  
    iatom += created_atoms;
    
  }


  if (iatom != totalatoms) {
    std::cerr << "Incorrect number of atoms specified in " << fname << std::endl;
    std::cerr << "atoms in file = " << iatom << ", total atoms = " << totalatoms << std::endl;
    return FORMAT_ERROR;
  }

  return SUCCESS;
}
