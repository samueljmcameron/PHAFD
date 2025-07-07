#include <algorithm>

#include "utility.hpp"

#include <iostream>
#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "fix.hpp"
#include "compute.hpp"
#include "dump.hpp"
#include "integrate.hpp"

#include "fix_output_file.hpp"

using namespace PHAFD_NS;

FixOutputFile::FixOutputFile(PHAFD *phafd) : Fix(phafd) {
  averaging = true;
};


void FixOutputFile::init(const std::vector<std::string> &v_line)
{

  compute = nullptr;
  fix = nullptr;
  
  Fix::init(v_line);

  std::vector<std::string> new_v_line(v_line);

  filename = new_v_line.at(1);
  every = std::stoi(new_v_line.at(2));
  
  if (every <= 0 )
    throw std::invalid_argument("Every must be >= 0.");  
  
  std::string arrname = new_v_line.at(3);
  output_type = new_v_line.at(4);
  
  std::string prefix = arrname.substr(0,2);
  std::string id = arrname.substr(2);
  
  if (prefix == "c_") {
    
    int index = utility::find_index(id,Compute::NAMES);
    compute = computes.at(index).get();
    
    local0start = compute->local0start;
    localNz = compute->localNz;
    localNy = compute->localNy;
    localNx = compute->localNx;
    
    utility::type_of_output(output_type,id,compute);

    
    
  } else if (prefix == "f_") {
    
    int index = utility::find_index(id,Fix::NAMES);
    fix = fixes.at(index).get(); 
    
    local0start = fix->local0start;
    localNz = fix->localNz;
    localNy = fix->localNy;
    localNx = fix->localNx;
    
    utility::type_of_output(output_type,id,compute);  

    
  } else
    throw std::runtime_error("Invalid input ID " + arrname +
			     std::string(" in compute/ft/spherical/bin."));



    
  
}


void FixOutputFile::setup()
{


  std::string tmp_name;

  if (compute != nullptr) {
    input_nc = compute->numberofcomponents;
    input_array_size = compute->array.size();
    input_array = compute->array.data();
    tmp_name =  compute->name;

  } else if (fix != nullptr) {
    input_nc = fix->numberofcomponents;
    input_array_size = fix->array.size();
    input_array = fix->array.data();
    tmp_name = fix->name;
  }

  if (output_type == "vector") {
    array.resize(input_array_size);
  }
  for (int i  = 0 ; i < array.size(); i++)
    array[i] = 0.0;
  

  if (commbrick->me == 0) {
    auto myfile = std::ofstream(filename);
    if (not myfile.is_open()) 
      throw std::runtime_error(std::string("Cannot open file ")
			       + filename);
    myfile << "# Output for " << name <<std::endl;

    myfile << "# TimeStep Number-of-rows" << std::endl;
    
    
    myfile << "# Row ";
    for (int col = 0; col < input_nc; col++)
      myfile << tmp_name << "[" << col << "] ";
    myfile << std::endl;
  
  }
}


void FixOutputFile::start_of_step()
{

  if (integrate->timestep % every == 0) {
    this_step = true;
    if (compute != nullptr) {
      compute->this_step = true;
      compute->start_of_step();
    } else if (fix != nullptr) {
      fix->this_step = true;
      fix->start_of_step();
    }
  }
  
  
}


void FixOutputFile::end_of_step() {


  if (integrate->timestep % every != 0) return;
  for (int i  = 0 ; i < array.size(); i++) {
    array[i] = 0.0;
  }

  if (compute != nullptr)	 {
    
    for (int i  = 0 ; i < array.size(); i++) {
      array[i] = input_array[i];
      
    }
    
  } else if (fix != nullptr) {
    
    for (int i  = 0 ; i < array.size(); i++) {
      array[i] = input_array[i];
      
    }
  }


  if (output_type == "vector") {

    if (commbrick->me == 0) {
      auto myfile = std::ofstream(filename,std::ios_base::app);
      if (not myfile.is_open()) 
	throw std::runtime_error(std::string("Cannot open file ")
				 + filename);
      
      int out_arr_rows = input_array_size/input_nc;      
      myfile << std::to_string(integrate->timestep) << " "
	     << std::to_string(out_arr_rows) << std::endl;


      for (int row = 0; row < out_arr_rows; row++) {
	myfile << row+1 << " ";
	for (int col = 0; col < input_nc; col++)
	  myfile << std::to_string(array[row+col*out_arr_rows]) << " ";
	myfile << std::endl;
      }
      
    }
  }


  if (compute != nullptr)	
    compute->this_step = false;
  else if (fix != nullptr)
    fix->this_step = false;

  this_step = false;

  
  /*
  int global_size;

  std::cout << "local_size = " << input_array_size << std::endl;
  MPI_Allreduce(&input_array_size,&global_size,1,MPI_INT,MPI_SUM,world);
  std::cout << "global_size = " << global_size << std::endl;  
  std::vector<double> output_array;
  //if (commbrick->me == 0) {
  output.resize(global_size);
    //std::cout << "global output size = " << output.size() << std::endl;
    //}



  // fill output_array only on processor zero
  MPI_Gather(input_array,input_array_size,MPI_DOUBLE,output.data(),
	     input_array_size,MPI_DOUBLE,0,world);

  if (commbrick->me == 0) {
    auto myfile = std::ofstream(filename,std::ios_base::app);
    if (not myfile.is_open()) 
      throw std::runtime_error(std::string("Cannot open file ")
			       + filename);
    
    myfile << std::to_string(integrate->timestep) << " "
	   << std::to_string(global_size) << std::endl;
    
    int count = 0;
    for (int i = 0; i < global_size; i++) {
      for (int dim = 0; dim < input_nc; dim++)
	myfile << std::to_string(output_array[count++]);
      myfile << std::endl;
      
    }

  }
  */
  
  
  
  

}
