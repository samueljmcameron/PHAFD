#include <iostream>
#include <algorithm>
#include <set>

#include "atom.hpp"
#include "utility.hpp"
#include "read_dump.hpp"
#include "read_atoms.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "group.hpp"

using namespace PHAFD_NS;

Atom::Atom(PHAFD *phafd)
  : Pointers(phafd),size_forward(3), size_reverse(3),size_border(6) {
  ntypes = -1;
  atomset = false;

};



//void Atom::init() {};

void Atom::setup(const std::vector<std::string> & v_line) {

  molecule_flag = sphere_flag = false;

  std::vector<int> Natomsperproc(commbrick->nprocs,0);

  
  int nargs = v_line.size();  
  int iarg = 0;

  while (v_line.at(iarg) != "natoms") {
    if (v_line.at(iarg) == "polymer") {
      molecule_flag = true;

    } else if (v_line.at(iarg) == "sphere") {
      sphere_flag = true;
    } else
      throw std::runtime_error("Unrecognised atom style.");


    iarg += 1;
    
  }

  iarg += 1;

  std::string tmpstr;

  int proc = 0;

  // get number of atoms per processor
  while (iarg < nargs) {

    tmpstr = v_line.at(iarg);
    std::string::size_type pos = tmpstr.find("*");
    if (pos != std::string::npos) {
      int atoms_per_proc = std::stoi(tmpstr.substr(0,pos));
      int nps = std::stoi(tmpstr.substr(pos+1));

      if (proc+nps > commbrick->nprocs)
	throw std::runtime_error("More atom groupings than there are processors.");
      for (int i = 0; i < nps; i++) {
	Natomsperproc.at(proc++) = atoms_per_proc;

      }
    } else
      Natomsperproc.at(proc++) = std::stoi(tmpstr);
    iarg += 1;
  }


  int Natoms = Natomsperproc.at(commbrick->me);


  

  nowned = 0;
  ngathered = 0;
  nlocal = 0;
  nghost = 0;

  nmols = 0;
  
  xs.resize(3,Natoms); 
  uxs.resize(3,Natoms); 
  Fs.resize(3,Natoms);
  tags.resize(Natoms);
  types.resize(Natoms);
  labels.resize(Natoms);
  images.resize(Natoms);

  if (molecule_flag) {
    moltags.resize(Natoms);
  }

  if (sphere_flag) {
    radius.resize(Natoms);
  }

  atomset = true;
  
}

void Atom::populate(const std::vector<std::string> & v_line) {

  if (v_line.at(0) == "read_atoms") {
     ReadAtoms read_atoms(phafd);
      
     int errflag = read_atoms.read_file(v_line.at(1));
     if (errflag != ReadAtoms::SUCCESS) errflag = 1;
     else errflag = 0;
     int total_errflag;
     MPI_Allreduce(&errflag,&total_errflag,1,MPI_INT,MPI_SUM,world);
     if (total_errflag)
       throw std::runtime_error("Could not read data atom file. ");

  } else if (v_line.at(0) == "read_dump") {

    ReadDump read_dump(phafd);

    std::vector<std::string> new_v_line(v_line);

    new_v_line.at(0) = "atom";
    read_dump.init(new_v_line);
    read_dump.process_attributes();

  } else
    throw std::runtime_error("Invalid atom_populate command.");

}

void Atom::check_tags_and_types()
/*
  Call this function after all atoms are created to check that things are as expected.

*/
{

  
  std::set<int> mset (moltags.begin(),moltags.end());
  std::set<int> tset(tags.begin(),tags.end());

  std::vector<int> mtmp (mset.begin(),mset.end());


  int errflag = 0;
  int totalerr;
  
  if (tags.size() != nowned) {
    std::cerr << "number of atom allocated on processor " << commbrick->me
	      << " doesn't match number of atoms created." << std::endl;
    errflag = 1;
  }

  MPI_Allreduce(&errflag,&totalerr,1,MPI_INT,MPI_MAX,world);
  if (totalerr)
    throw std::runtime_error("Incorrect atoms allocated on processor(s).");



  if (tset.size() != tags.size()) {
    std::cerr << "duplicate IDs on processor " << commbrick->me << "." << std::endl;
    errflag = 1;
  }

  MPI_Allreduce(&errflag,&totalerr,1,MPI_INT,MPI_MAX,world);
  if (totalerr)
    throw std::runtime_error("duplicate IDs found on processor(s).");



  // check if any IDs are the same across all processors
  utility::check_MPI_duplicates(tags,world,commbrick->me,commbrick->nprocs,"IDs");

  // check if any molIDs are the same across processors
  utility::check_MPI_duplicates(mtmp,world,commbrick->me,commbrick->nprocs,"molIDs");



  
  
  int typemax;

  if (types.size() == 0)
    typemax = -1;
  else
    typemax = *(std::max_element(types.begin(),types.end()));



  int all_typemax; 

  MPI_Allreduce(&typemax,&all_typemax,1,MPI_INT,MPI_MAX,world);

  ntypes = all_typemax + 1;

  groups.push_back(std::make_unique<Group>(phafd));
  groups.back()->create_all();

  if (molecule_flag) {
    std::set<int> molset(moltags.begin(),moltags.end());
    nmols = molset.size();

  }

  

}


int Atom::add_atom(std::vector<std::string> v_line, int iatom)
/*
  Create a single sphere at index iatom. Input should have format id type x 

 */
{
  

  int l_index = 0;

  tags[iatom] = std::stoi(v_line.at(l_index++));

  if (molecule_flag) {
    moltags[iatom] = std::stoi(v_line.at(l_index++));
  }

  
  types[iatom] = std::stoi(v_line.at(l_index++));
  xs(0,iatom) = std::stod(v_line.at(l_index++));
  xs(1,iatom) = std::stod(v_line.at(l_index++));
  xs(2,iatom) = std::stod(v_line.at(l_index++));

  if (sphere_flag)
    radius[iatom] = std::stod(v_line.at(l_index++));
  
  labels[iatom] = Atom::OWNED;
  images[iatom] = domain->set_image();

  nowned += 1;
  
  return 1;
  

}


int Atom::unpack_reverse(int n,const std::vector<int> & sendlist, double *buf)
{

  int m = 0;
  int j;
  for (int i = 0; i < n; i++) {
    j = sendlist[i];
    Fs(0,j) += buf[m++];
    Fs(1,j) += buf[m++];
    Fs(2,j) += buf[m++];
  }

  return m;
  
}

int Atom::pack_comm(int n,const std::vector<int> & sendlist, double *buf,
		    const std::array<double,3> & pbc)
{

  int m = 0;
  int j;
  for (int i = 0; i < n; i++) {
    j = sendlist[i];
    buf[m++] = xs(0,j) + pbc[0];
    buf[m++] = xs(1,j) + pbc[1];
    buf[m++] = xs(2,j) + pbc[2];
  }

  return m;
  
}

int Atom::pack_comm_in_z(int n,const std::vector<int> & sendlist, double * buf,
			 const std::vector<double> &pbcz )
{
 int m = 0;
  int j;
  for (int i = 0; i < n; i++) {
    j = sendlist[i];
    buf[m++] = xs(0,j);
    buf[m++] = xs(1,j);
    buf[m++] = xs(2,j) + pbcz[i];
  }
  return m;

}

int Atom::pack_border_in_z(int nsend, const std::vector<int> & sendlist,double *buf,
			   const std::vector<double> & pbcz,
			   const std::vector<int> & labelz)
{

  int j;

  int m = 0;

  for (int i = 0; i < nsend; i++) {
    j = sendlist[i];
    buf[m++] = xs(0,j);
    buf[m++] = xs(1,j);
    buf[m++] = xs(2,j) + pbcz[i];
    buf[m++] = ubuf(tags[j]).d;
    buf[m++] = ubuf(types[j]).d;
    buf[m++] = ubuf(labelz[i]).d;

  }

  return m;
}


int Atom::pack_border(int nsend, const std::vector<int> & sendlist,double *buf,
		      const std::array<double,3> & pbc)
{

  int j;

  int m = 0;

  for (int i = 0; i < nsend; i++) {
    j = sendlist[i];
    buf[m++] = xs(0,j) + pbc[0];
    buf[m++] = xs(1,j) + pbc[1];
    buf[m++] = xs(2,j) + pbc[2];
    buf[m++] = ubuf(tags[j]).d;
    buf[m++] = ubuf(types[j]).d;
    buf[m++] = ubuf(Atom::GHOST).d;
    
  }

  return m;
}



void Atom::unpack_border(int nrecv,int first, const double *buf)
{
  int m = 0;
  int last = first + nrecv;

  xs.conservativeResize(Eigen::NoChange,last);
  Fs.conservativeResize(Eigen::NoChange,last);
  tags.resize(last);
  types.resize(last);
  labels.resize(last);

  
  for (int i = first; i < last; i++) {

    xs(0,i) = buf[m++];
    xs(1,i) = buf[m++];
    xs(2,i) = buf[m++];
    tags[i] = (int) ubuf(buf[m++]).i;
    types[i] = (int) ubuf(buf[m++]).i;
    labels[i] = (int) ubuf(buf[m++]).i;
    
  }

  
  return;
}
