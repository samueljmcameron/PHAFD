
#include "utility.hpp"
#include "phafd.hpp"
#include "fftw_arr/array3d.hpp"
#include "compute.hpp"
#include "fix.hpp"
#include <random>

#include <algorithm>
#include <exception>
#include <set>

#include <iostream>

std::vector<std::string> PHAFD_NS::utility::split_line(std::string& line)
/*
  Given an input string, split it into words (separated by whitespace) and
  return the vector of these words - removing any training comments (starting with
  '#').

 */
{
  // remove any trailing comments
  line = line.substr(0,line.find("#",0));
  
  std::vector<std::string> out;
  
  std::string tmp;
  
  
  std::size_t index;
  
  index = 0;
  while(index < line.size()) {
    
    for (; index != line.size()
	   && !isspace(line[index]); ++index)
      tmp += line[index];
    
    if (tmp != "") {
      out.push_back(tmp);
    }
    tmp = "";
    
    index += 1;
  }
  
  return out;
}


void PHAFD_NS::utility::replacePercentages(std::string &raw,int id)
/*
  Given a string which has "%" in its name, replace this with the
  the integer id.   E.g. if id = 4 and raw is

  "fname_p%.txt"

  the output would be

  "fname_p4.txt"
  
 */
{

  std::string::size_type vstart;

  while (true) {
    
    vstart = raw.find("%");
    
    if (vstart != std::string::npos) 
      raw.replace(vstart,1,std::to_string(id));
    else break;
  }
  
  return;
}
  
void PHAFD_NS::utility::convertVariables(std::string &raw,
					 std::map<std::string, std::string> const& varMap)
/*

  Given a string, replace any set of characters with the form ${var} to the
  value of var, where var and its value must be in the map varMap. E.g. for
  a var map {"name" : "Jim"}, the string

  "Hello my name is ${name}othy - well actually it's just ${name}."
  
  would convert to
  
  "Hello my name is Jimothy - well actually it's just Jim."

  
 */
{
  
  
  
  std::string::size_type vstart,vend;
  
  vend = 0;
  while (vend != raw.size()) {
    
    vstart = raw.find("${");
    
    if (vstart != std::string::npos) {
      
      vend = raw.find("}",vstart+2);
      
      if (vend != std::string::npos) {
	std::string tmp = raw.substr(vstart+2,vend-vstart-2);
	if (tmp == "")
	  throw std::runtime_error("No content ('${}') in input file.");
	bool found_key = false;
	
	for (const auto &xm: varMap) {
	  if (xm.first == tmp) {
	    found_key = true;
	    raw.erase(vstart,vend - vstart+1);
	    raw.insert(vstart,xm.second);
	  }
	}
	if (!found_key) {
	  std::string errorMessage
	    = std::string("No such command line variable ") + tmp;
	  throw std::runtime_error(errorMessage);
	}
	
	
	
      } else {
	throw::std::runtime_error("Missing '}' in input script.");
      }
      
    } else {
      vend = raw.size();
    }
  }
  
  return;
}


void PHAFD_NS::utility::check_MPI_duplicates(const std::vector<int> & vec,
					     MPI_Comm comm,int id,int mpi_size,
					     std::string vecname)
/*
  check that atom IDs are not duplicated across processors, and that all vectors are filled.
  input vector should not have duplicate values (on same processor) but does not need to be
  sorted.
*/ 
{

  
  std::vector<int> number_to_check;

  if (id == 0)
    number_to_check.resize(mpi_size);

  
  int sendsize = vec.size();
  
  // gather all sendsizes to process 0
  MPI_Gather(&sendsize,1,MPI_INT,number_to_check.data(),1,MPI_INT,0,comm);
  // so now process 0 has number_to_check = [x,y,z,...] where x is size of vector on p0,
  //  y is size of vector on p1, ... etc.

  
  std::vector<int> vec_to_check; // vector to be checked for duplicates (only filled on process 0)
  std::vector<int> displs; // vector which specifies the offsets required for MPI_Gatherv


  
  if (id == 0) {
    
    displs.resize(mpi_size);
    
    int totalsize = 0;
    int count = 0;
    for (int num : number_to_check) {
      displs[count++] = totalsize;
      totalsize += num;
    }

    vec_to_check.resize(totalsize);
  }



  MPI_Gatherv(&vec[0],sendsize,MPI_INT,vec_to_check.data(),number_to_check.data(),
	      displs.data(),MPI_INT,0,comm);  

  int intersections = 0;
  
  if (id == 0) {


    std::set<int> vec_set(vec_to_check.begin(),vec_to_check.end());

    intersections = vec_to_check.size()-vec_set.size();
  }

  MPI_Bcast(&intersections,1,MPI_INT,0,comm);


  
  if (intersections) throw std::runtime_error(std::to_string(intersections)
					      + std::string(" duplicate ") + vecname
					      + std::string(" across processors."));
  
  return;  
}


int PHAFD_NS::utility::make_unique_seed(int baseseed,const MPI_Comm &comm,
					int id, int nprocs)
{

  std::mt19937 gen;
  std::uniform_int_distribution<int> integer_dist;

  std::vector<int> processor_seeds(nprocs);

  
  gen.seed(baseseed);

  if (id == 0) {
    for (auto & num : processor_seeds)
      num = integer_dist(gen);
  }

  MPI_Bcast(processor_seeds.data(),nprocs,MPI_INT,0,comm);


  return processor_seeds.at(id);

}



int PHAFD_NS::utility::find_brackets(std::string &id)
/*

  Check if "[%d]" is in string. If it is, return %d. If not,
  return -1. Also removes "[%d]" from the string (if it exists).


 */
{
  std::size_t pos = id.find("[");
      
  int arr_comp_num;
	
  if (pos != std::string::npos) {
    
    std::size_t spos = id.find("]");
    if (spos == std::string::npos)
      throw std::runtime_error("No closing bracket on ID "
			       + id);
    arr_comp_num = std::stoi(id.substr(pos,spos-pos));
    id = id.substr(0,pos);
  } else
    throw std::runtime_error("Expected bracketed expression in ID "
			     + id);


  return arr_comp_num;
}


int PHAFD_NS::utility::
find_index(std::string id,const std::vector<std::string> &names)
{
  
  int index = 0;
  
  for (auto &name : names) {
    
    if (id == name) {
      break;
    }
    index += 1;
  }
  if (index == names.size())
    throw std::runtime_error("ID " + id
			     + std::string("doesn't exist."));

  return index;
}


template <class T>
void PHAFD_NS::utility::type_of_output(std::string check_type,
				       std::string id,T * cls)
{
  bool anerror = false;
  if (check_type == "grid") {


    if (!cls->per_grid)
      anerror = true;
  
  } else if (check_type == "ftgrid") {
    if (!cls->per_ftgrid)
      anerror = true;
  
  } else if (check_type == "atom") {
    if (!cls->per_atom)
      anerror = true;
  
  } else if (check_type == "vector") {
    if (!cls->vector)
      anerror = true;
  
  } else if (check_type == "scalar") {
    if (!cls->scalar)
      anerror = true;


    
  } else if (check_type != "")
    throw std::runtime_error("invalid check_type on ID " + id);
  if (anerror)
    throw std::runtime_error("ID " + id +
			     std::string(" does not have ")
			     + check_type + std::string(" quantity."));
  
}

void PHAFD_NS::utility::
find_array_component(std::string arrname, PHAFD *phafd,
		     fftwArr::array3D<double> * array,
		     std::string check_output)
{
  int arr_comp_num = find_brackets(arrname);
  
  if (arrname.rfind("c_",0) == 0) {

    std::string cid = arrname.substr(2);

    int index = find_index(std::string(cid),Compute::NAMES);

    auto cmp = phafd->computes.at(index).get();

    PHAFD_NS::utility::type_of_output(check_output,cid,cmp);
    array = cmp->realFFTWarray.at(arr_comp_num);
    
  } else if (arrname.rfind("f_",0) == 0) {
    
    std::string fid = arrname.substr(2);
    
    int arr_comp_num = find_brackets(fid);

    int index = find_index(std::string(fid),Fix::NAMES);

    auto fx = phafd->fixes.at(index).get();
    
    PHAFD_NS::utility::type_of_output(check_output,fid,fx);
    array = fx->realFFTWarray.at(arr_comp_num);

    
  } else
    array = nullptr;

  return;
}


void PHAFD_NS::utility::
find_array_component(std::string arrname, PHAFD *phafd,
		     fftwArr::array3D<std::complex<double>> * array,
		     std::string check_output)
{

  int arr_comp_num = find_brackets(arrname);
  
  if (arrname.rfind("c_",0) == 0) {
    
    std::string cid = arrname.substr(2);
    
    int index = find_index(std::string(cid),Compute::NAMES);

    auto cmp = phafd->computes.at(index).get();
    
    PHAFD_NS::utility::type_of_output(check_output,cid,cmp);    
    array = cmp->complexFFTWarray.at(arr_comp_num);
    
    
  } else if (arrname.rfind("f_",0) == 0) {
    
    std::string fid = arrname.substr(2);
      
    int index = find_index(std::string(fid),Fix::NAMES);
    
    auto fx = phafd->fixes.at(index).get();    
    
    PHAFD_NS::utility::type_of_output(check_output,fid,fx);
    array = fx->complexFFTWarray.at(arr_comp_num);
    
    
  } else
    array = nullptr;
  
  return;
}

template void PHAFD_NS::utility::type_of_output<
  PHAFD_NS::Compute>(std::string ,std::string ,Compute *);

template void PHAFD_NS::utility::type_of_output<
  PHAFD_NS::Fix>(std::string , std::string ,Fix *);
//template void PHAFD_NS::utility::find_array_component<
//  fftwArr::array3D<double>>(const std::string &,PHAFD *,
//			    fftwArr::array3D<double> *);

//template void PHAFD_NS::utility::find_array_component<
//  fftwArr::array3D<std::complex<double>
//		   >>(const std::string &,PHAFD *,
//		      fftwArr::array3D<std::complex<double>> *);
