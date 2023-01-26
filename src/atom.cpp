#include <iostream>
#include <algorithm>
#include <set>

#include "atom.hpp"
#include "utility.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"

#include "beadrodpmer/iovtk.hpp"

using namespace PHAFD_NS;

Atom::Atom(PHAFD *phafd)
  : Pointers(phafd),size_forward(3), size_reverse(3),size_border(6) {};



//void Atom::init() {};

void Atom::setup(int Natoms,std::vector<std::string> styles) {

  molecule_flag = sphere_flag = false;

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

  if (std::find(styles.begin(), styles.end(), "polymer")!=styles.end()) {
    molecule_flag = true;
    moltags.resize(Natoms);
  }

  if (std::find(styles.begin(), styles.end(), "sphere")!=styles.end()) {
    sphere_flag = true;
    radius.resize(Natoms);
  }
  
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

  ntypes = 2;

}



int Atom::add_polymer(std::vector<std::string> v_line,int startatom)
/* 
   Input line should have the form "idstart-idend mol t1*l1 t2*l2 ... tn*ln" where
   idstart is the id of the first atom in the polymer, idend-1 is the id of the last atom
   in the polymer (so y-x is the polymer length), mol is the molecule id,
   t1 is the type of the first l1 atoms in the polymer, etc,
   so l1 + l2 + ... + ln = idstart-idend .
*/
{

  if (!molecule_flag)
    throw std::invalid_argument("Must have molecular style atoms for creating polymer.");    
  
  if (v_line.size() < 3)
    throw std::invalid_argument("Invalid arguments for creating polymer.");    

  std::size_t pos = v_line.at(0).find("-");

  if (pos == std::string::npos)
    throw std::invalid_argument("Invalid ID specification when creating polymer.");


  int idstart,idend,molid;


  idstart = std::stoi(v_line.at(0).substr(0,pos));

  idend = std::stoi(v_line.at(0).substr(pos+1));

  molid = std::stoi(v_line.at(1));

  
  v_line.erase(v_line.begin(),v_line.begin()+2);

  std::vector<int> spectypes,lengthtypes;

  std::string fname = v_line.back();

  v_line.pop_back();
  
  for (auto word : v_line) {
    pos = word.find("*");
    if (pos == std::string::npos)
      throw std::invalid_argument("Invalid type specification when creating polymer.");

    spectypes.push_back(std::stoi(word.substr(0,pos)));

    lengthtypes.push_back(std::stoi(word.substr(pos+1)));

  }


  int lt_sum = 0;

  for (int lt : lengthtypes) lt_sum += lt;

  if (lt_sum != idend - idstart)
    throw std::invalid_argument("Type specification does not match number of "
				"atoms in creating polymer.");

  int lt_index = 0;
  lt_sum = lengthtypes.at(0);



  // set IDs, molIDs, types
  for (int iatom = startatom; iatom < idend-idstart+startatom; iatom ++) {

    if (iatom >= tags.size()) {
      throw std::invalid_argument("atom index " + std::to_string(iatom)
				  + std::string(" out of range when creating polymer")
				  + std::string(" starting at index ")
				  + std::to_string(startatom));
    }
    tags[iatom] = idstart + iatom - startatom;


    
    moltags[iatom] = molid;

    if (iatom - startatom == lt_sum) {
      lt_sum += lengthtypes[++lt_index];
    }
    
    types[iatom] = spectypes[lt_index];
    images[iatom] = domain->set_image();
    labels[iatom] = Atom::OWNED;
    
  }


  nowned += idend - idstart;

  nmols += 1;

  BeadRodPmer::ioVTK::readVTKPolyData(xs.middleCols(startatom,idend-idstart),fname);

  return idend - idstart;
}


int Atom::add_sphere(std::vector<std::string> v_line, int iatom)
/*
  Create a single sphere at index iatom. Input should have format id type x 

 */
{
  

  int l_index = 0;

  tags[iatom] = std::stoi(v_line.at(l_index++));

  if (molecule_flag) {
    moltags[iatom] = std::stoi(v_line.at(l_index++));
    nmols += 1;
  }

  
  types[iatom] = std::stoi(v_line.at(l_index++));
  xs(0,iatom) = std::stod(v_line.at(l_index++));
  xs(1,iatom) = std::stod(v_line.at(l_index++));
  xs(2,iatom) = std::stod(v_line.at(l_index++));
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
