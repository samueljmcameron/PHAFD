
#include "utility.hpp"

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
