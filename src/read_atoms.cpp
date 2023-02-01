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


  
  // set per-style vectors
  errflag =  read_perStyle();


  if (errflag != SUCCESS) errflag = 1;
  else errflag = 0;

  MPI_Allreduce(&errflag,&total_errflag,1,MPI_INT,MPI_SUM,world);
  
  if (total_errflag) return FORMAT_ERROR;



  // allocate the space needed in Atom
  
  int totalatoms = 0;

  for (auto num : atomsPerStyle)
    totalatoms += num;

    
  atoms->setup(totalatoms,styles);


  // create the atoms from the file info
  
  errflag = create_atoms();
  
  MPI_Allreduce(&errflag,&total_errflag,1,MPI_INT,MPI_SUM,world);

  if (total_errflag) return FORMAT_ERROR;
  
  datafile.close();
  datafile.clear();
  
  return SUCCESS;
  
}





int ReadAtoms::read_perStyle()
/*
  Read in per-style vectors styles, particlesPerStyle, and atomsPerStyle.
*/
{

  int errflag;
  
  // set styles from "particle_style" command
  errflag = reset_styles(); 
  
  if (errflag == FORMAT_ERROR) 
    return errflag;
  
  // set # of particles per style from "particlesperstyle" command
  errflag = reset_particlesPerStyle(); 

  if (errflag == FORMAT_ERROR)
    return errflag;
  
  
  // set # of atoms per style from particle info and "atomsperpolymer" command
  errflag = reset_atomsPerStyle(); 
  
  if (errflag == FORMAT_ERROR)

    return errflag;

  return SUCCESS;
}




int ReadAtoms::create_atoms()
/* Construct atoms from the lines of the data file. */
{

  // allocate unique_ptr for atoms



  int iatom = 0;


  std::string line;

  std::vector<std::string> v_line;

  int created_atoms;
  
  for (int istyle = 0; istyle < styles.size(); istyle++) {

    int imol = 0;

    while (std::getline(datafile,line)) {
      
      if (line == "" || line[0] == '#') continue;

      v_line = utility::split_line(line);
    
      if (styles[istyle] == "polymer") {

	try {
	  
	  created_atoms = atoms->add_polymer(v_line,iatom);
	  
	} catch (std::invalid_argument &inv) {
	  
	  std::cerr << inv.what();
	  return FORMAT_ERROR;
	  
	}
	
      } else if (styles[istyle] == "sphere") {

	try {
	  
	  created_atoms = atoms->add_sphere(v_line,iatom);
	  
	} catch (std::invalid_argument &inv) {
	  
	  std::cerr << inv.what();
	  return FORMAT_ERROR;
	  
	}

      } else {
	std::cerr << "Invalid particle style." << std::endl;
	return FORMAT_ERROR;
      }

      imol += 1;
      iatom += created_atoms;
      
      if (imol ==  particlesPerStyle[istyle]) break;





    }

    
  }


  atoms->check_tags_and_types();

  return SUCCESS;
}



int ReadAtoms::reset_styles() {
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


int ReadAtoms::reset_particlesPerStyle()
/*
  (Re)set the number of particles on the processor. Must be called after reset_styles()
  method.
*/
{

  particlesPerStyle.clear();

  std::string line;

  
  std::vector<std::string> v_line; // to store words from line of file

  
  while (std::getline(datafile,line)) {
    
    if (line == "" || line[0] == '#') continue;
    
    
    break;
  }

  v_line = utility::split_line(line);
  
  if (v_line.size() != 1 + styles.size()
      || v_line.at(0) != "particlesperstyle" ) {
    std::cerr << "must specify 'particlesperstyle' in " << fname << std::endl;
    return FORMAT_ERROR;
  }
  
  
  v_line.erase(v_line.begin());

  for (auto elem : v_line) {
    
    try {
      particlesPerStyle.push_back(std::stoi(elem));
    } catch (std::invalid_argument &err) {
      std::cerr << "'particlesperstyle' requires integer values in " << fname << std::endl;
      return FORMAT_ERROR;
    }
    
  }
  
  return SUCCESS;
}


int ReadAtoms::reset_atomsPerStyle()
/*
  (Re)set the number of atoms on the processor. Must be called after
  reset_particlesPerStyle() method.
*/
{

  atomsPerStyle.resize(styles.size());
  
  for (auto istyle = 0; istyle < styles.size(); istyle++) {
    
    atomsPerStyle[istyle] = 0;
    
    if (styles[istyle] == "polymer") {

      atomsPerPolymer.clear();
      
      std::string line;
      std::vector<std::string> v_line;
      
      while (std::getline(datafile,line)) {
	
	if (line == "" || line[0] == '#') continue;
	
	break;
      }
      
      v_line = utility::split_line(line);
      
      if (v_line.at(0) != "atomsperpolymer") {
	std::cerr << "'polymer' style requires atomsperpolymer to be specified in "
		  << fname << std::endl;	
	return FORMAT_ERROR;
      }
      
      v_line.erase(v_line.begin());
      
      if (v_line.size() != particlesPerStyle[istyle]) {
	std::cerr << "'atomsperpolymer' must be specified for each polymer in "
		  << fname << std::endl;	
	return FORMAT_ERROR;
      }

      for (auto snum : v_line) {
	try {

	  atomsPerPolymer.push_back(std::stoi(snum));
	} catch (std::invalid_argument &inv) {
	  std::cerr << "'atomsperpolymer' requires integer values in "
		    << fname << std::endl;	
	  return FORMAT_ERROR;
	}
	atomsPerStyle[istyle] += atomsPerPolymer.back();
      }

      
    } else {
      atomsPerStyle[istyle] = particlesPerStyle[istyle];
    }
    
  }

  return SUCCESS;
}
