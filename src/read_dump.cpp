#include <iomanip>
#include <iostream>
#include "read_dump.hpp"

#include "domain.hpp"
#include "grid.hpp"
#include "fix.hpp"
#include "compute.hpp"
#include "atom.hpp"
#include "fftw_arr/array3d.hpp"
#include "integrate.hpp"

using namespace PHAFD_NS;

/* dump file structure is going to be that each timestep is written to an individual (VTK)
   file, and then there is a .pvd file to display the timesteps collectively. Files are
   always going to be outputted in parallel. */
ReadDump::ReadDump(PHAFD *phafd) : Pointers(phafd) {
};



void ReadDump::init(const std::vector<std::string> &v_line) {

  
  dump_type = v_line.at(0);


  filename = v_line.at(1);

  auto it = v_line.begin()+2;
  while (it != v_line.end()) {
    
    if (*it == "attributes") {
      ++it;
      int numattributes = std::stoi(*it);
      ++it;
      for (int i = 0; i < numattributes; i++) {
	if (it == v_line.end())
	  throw std::invalid_argument("Incorrect number of attributes in dump file.");
	attributes.insert(*it);
	++it;
      }

      if (attributes.size() != numattributes)
	throw std::invalid_argument("Duplicate attribute values in dump file.");
    } else
      throw std::invalid_argument("Third word in dump must be 'attributes'.");
  }

  if (attributes.size() == 0)
    throw std::invalid_argument("Need attributes to read dump files.");

  for (auto &word : attributes) {
    if (word.rfind("c_",0) == 0 || word.rfind("f_",0) == 0) {

      throw std::runtime_error("Cannot read computes from dump");


      
    } 
  }
  
}



void ReadDump::process_attributes()
{
  if (dump_type == "atom") {
    process_atom_attributes();
  } else if (dump_type == "grid") {
    process_grid_attributes();
  } else {
    throw std::runtime_error("Can only read atoms or grid from dump file.");
  }
}

void ReadDump::process_grid_attributes()
/* Start from the beginning of the dump file. */
{
  int localerr = 0;
  int globalerr;
  auto myfile = std::fstream(filename, std::ios::in | std::ios::binary);
  if (myfile.fail()) {
    std::cerr << "could not open file " << filename << std::endl;
    localerr = 1;
  }
  MPI_Allreduce(&localerr,&globalerr,1,MPI_INT,MPI_SUM,world);

  if (globalerr)
    throw std::runtime_error("file not found.");
  std::string line,expected;


  // necessary attributes
  std::vector<std::string> necessary = {"phi"};

  for (auto &word : necessary)
    if (attributes.find(word) == attributes.end())
      throw std::runtime_error("Read dump needs attribute " + word + ".");

  
  std::map<std::string,int> offsets;
  while (std::getline(myfile,line)) {

    if (line == "<AppendedData encoding=\"raw\">")
      break;

    for (auto &word : attributes) {
      expected = "<DataArray Name=\"" + word + "\" type=\"Float64\" format=\"appended\" offset=\"";
      if (line.substr(0,expected.length()) == expected) {
	offsets.insert({word,std::stoi(line.substr(expected.length()))});
	break;
      }
    }
  }

  // read in the underscore
  char memblock[2];
  myfile.read(memblock,1);

  // then read in all the array files


  for (const auto& [key,value] : offsets) {
    if (key == "phi") {
      read_binary_data(myfile,grid->phi.get());
    } else if (key == "chempot") {
      read_binary_data(myfile,grid->chempot.get());
    } else if (key == "gradphi_x") {
      read_binary_data(myfile,grid->gradphi[0].get());
    } else if (key == "gradphi_y") {
      read_binary_data(myfile,grid->gradphi[1].get());
    } else if (key == "gradphi_z") {
      read_binary_data(myfile,grid->gradphi[2].get());
    } else {
      throw std::runtime_error("ReadDump error: Attribute does not exist.");
    }
  }
  return;
}


void ReadDump::process_atom_attributes()
/* Start from the beginning of the dump file. */
{

  int localerr = 0;
  int globalerr;
  auto myfile = std::fstream(filename, std::ios::in);
  if (myfile.fail()) {
    std::cerr << "could not open file " << filename << std::endl;
    localerr = 1;
  }
  MPI_Allreduce(&localerr,&globalerr,1,MPI_INT,MPI_SUM,world);

  if (globalerr)
    throw std::runtime_error("file not found.");

  std::string line;
  std::set<std::string> attributes_copy;

  // necessary attributes
  std::vector<std::string> necessary = {"x","type", "ID"};

  if (atoms->molecule_flag) {
    necessary.push_back("molID");
    necessary.push_back("ix");
  }
  if (atoms->sphere_flag)
    throw std::runtime_error("Not yet implemented read dump for radius.");

  for (auto &word : necessary)
    if (attributes.find(word) == attributes.end())
      throw std::runtime_error("Read dump needs attribute " + word + ".");
  
  
  while (std::getline(myfile,line)) {

    for (auto &word : attributes) {
      if (line ==  "<DataArray Name=\"" + word + "\" type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">") {
	std::getline(myfile,line);

	if (word == "x")
	  read_ascii_data(line,atoms->xs);
	else if (word == "F")
	  read_ascii_data(line,atoms->Fs);
	else if (word == "ux")
	  read_ascii_data(line,atoms->uxs);
	else
	  throw std::runtime_error("Attribute " + word + " not found in file.");
	attributes_copy.insert(word);
	break;
      }	else if (line ==  "<DataArray Name=\"" + word + "\" type=\"Int64\" NumberOfComponents=\"1\" format=\"ascii\">") {
	std::getline(myfile,line);
	
	if (word == "ix")
	  read_ascii_data(line,atoms->images);
	else if (word == "ID") 
	  read_ascii_data(line,atoms->tags);
	else if (word == "molID")
	  read_ascii_data(line,atoms->moltags);
	else if (word == "type") 
	  read_ascii_data(line,atoms->types);
	else if (word == "label")
	  read_ascii_data(line,atoms->labels);
	else
	  throw std::runtime_error("Attribute " + word + " not found in file.");
	attributes_copy.insert(word);
	
	break;
      } 
      
      
    }
  }

  atoms->nowned = atoms->tags.size();


  if (attributes_copy != attributes) {
    std::string errstring = "Only found attributes ";
    for (auto &word : attributes_copy)
      errstring += word  + ", ";
    throw std::runtime_error(errstring);
  }
  
  return;
}
 


void ReadDump::read_ascii_data(const std::string &line,
			       Eigen::Ref<Eigen::Matrix3Xd> array)
{


  std::stringstream ss(line);
  for (int i = 0; i < array.cols(); i++) {
    ss >> array.col(i)(0);
    ss >> array.col(i)(1);
    ss >> array.col(i)(2);
    
  }

}


void ReadDump::read_ascii_data(const std::string &line,
			       std::vector<int> & array)
{

  std::stringstream ss(line);
  for (int i = 0; i < array.size(); i++) 
    ss >> array[i];

}


void ReadDump::read_binary_data(std::fstream &myfile,
				fftwArr::array3D<double> *array) {

  unsigned int bytelength;

  myfile.read((char*)&bytelength,sizeof(bytelength));

  if (bytelength != array->Nz()*array->Ny()*array->Nx()*sizeof(double))
    throw std::runtime_error("incorrect size " + array->get_name() + " in file.");
  
  // since real fftw arrays aren't contiguous, need to read each row separately.
  
  for (int i = 0; i < array->Nz(); i++) {
    for (int j = 0; j < array->Ny(); j++) {
      myfile.read((char*)&(*array)(i,j,0),sizeof(double)*array->Nx());

    }
  }

  return;

}

