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

  // set styles

  std::vector<std::string> styles; // names of the particles styles (size 1 if not hybrid)

  errflag = reset_styles(styles);   



  if (errflag != SUCCESS) errflag = 1;
  else errflag = 0;

  MPI_Allreduce(&errflag,&total_errflag,1,MPI_INT,MPI_SUM,world);
  
  if (total_errflag) return FORMAT_ERROR;


  int totalatoms; // total number of atoms
  errflag = get_total_atoms(totalatoms);


  if (errflag != SUCCESS) errflag = 1;
  else errflag = 0;

  MPI_Allreduce(&errflag,&total_errflag,1,MPI_INT,MPI_SUM,world);

  if (total_errflag) return FORMAT_ERROR;

  
  atoms->setup(totalatoms,styles);


  // create the atoms from the file info
  
  errflag = create_atoms(totalatoms);
  
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
    return FORMAT_ERROR;
  }

  atoms->check_tags_and_types();
  std::cout << atoms->ntypes << std::endl;

  return SUCCESS;
}



int ReadAtoms::reset_styles(std::vector<std::string> & styles) {
  /*
    find the first non-empty, non-commented line of the file,
    which must start with "particle_style". Then store the particle style name(s)
    in the vector styles.
    
  */

  styles.clear();
  
  std::string line;
  
  std::vector<std::string> v_line; // to store words from line of file
  
  while (std::getline(datafile,line)) {
    
    if (line == "" || line[0] == '#') continue;
    
    v_line = utility::split_line(line);
    
    
    if (v_line.size() < 2 || v_line.at(0) != "particle_style") {
      std::cerr << "first line of " << fname << "must 'particle_style (hybrid) style(1) ...'"
		<< std::endl;
	  
      return FORMAT_ERROR;
    }
    
    if (v_line.at(1) == "hybrid") {
      if (v_line.size() < 3) {
	std::cerr << "hybrid style needs multiple substyles in " << fname 
		  << std::endl;
	return FORMAT_ERROR;
      }
      
      v_line.erase(v_line.begin(),v_line.begin()+2);
    } else {
      if (v_line.size() != 2) {
	std::cerr << "only one style permitted (if you need multiple styles, use hybrid keyword)"
		  << fname << std::endl;
	return FORMAT_ERROR;
      }
      
      v_line.erase(v_line.begin(),v_line.begin()+1);
    }
    
    std::set<std::string> tmpstyles;
    for (auto pstyle : v_line) {
      
      if (pstyle != "point" && pstyle != "sphere" && pstyle != "polymer") {
	std::cerr << "invalid style name in " << fname << std::endl;

	return FORMAT_ERROR;
      }
      
      styles.push_back(pstyle);
      tmpstyles.insert(pstyle);
    }
    
    if (tmpstyles.size() != styles.size()) {
      std::cerr << "Cannot have duplicate particle styles in " << fname << std::endl;
      return FORMAT_ERROR;
    }
    
    break;
  }

  return SUCCESS;
  
};


int ReadAtoms::get_total_atoms(int &totalatoms)
{

  std::string line;
  
  std::vector<std::string> v_line; // to store words from line of file

  
  totalatoms = 0;
  while (std::getline(datafile,line)) {
    
    if (line == "" || line[0] == '#') continue;
    
    v_line = utility::split_line(line);
    if (v_line.at(0) != "total_atoms") {
      std::cerr << "Expected 'total_atoms' command in " << fname << std::endl;
      return FORMAT_ERROR;
    }

    totalatoms = std::stoi(v_line.at(1));

    if (totalatoms < 0) {
      std::cerr << "'total_atoms' must be non-negative in " << fname << std::endl;
      return FORMAT_ERROR;
    }

    break;
  }

  return SUCCESS;
  
}
