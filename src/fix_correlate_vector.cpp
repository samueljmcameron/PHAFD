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

#include "fix_correlate_vector.hpp"

using namespace PHAFD_NS;

FixCorrelateVector::FixCorrelateVector(PHAFD *phafd) : Fix(phafd) {
  vector = true;
};


void FixCorrelateVector::init(const std::vector<std::string> &v_line,
			      bool add_to_names)
{


  compute = nullptr;
  fix = nullptr;
  
  Fix::init(v_line,add_to_names);

  std::vector<std::string> new_v_line(v_line);

  every = std::stoi(new_v_line.at(1));
  repeat = std::stoi(new_v_line.at(2));
  freq = std::stoi(new_v_line.at(3));
  std::string arrname = new_v_line.at(4);
  output_type = new_v_line.at(5);

  if (output_type != "vector")
    throw std::invalid_argument("Need to specify vector output type");


  if (every <= 0 || repeat <= 0 || freq <= 0)
    throw std::invalid_argument("Every, repeat, and freq must all be >= 0.");  
  if (every*(repeat-1) > freq)
    throw std::invalid_argument("Every*(repeat-1) exceeds freq in fixgridave.");
  if (freq % every != 0)
    throw std::invalid_argument("Every and repeat exceed freq in fixgridave.");  

  
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
    
    utility::type_of_output(output_type,id,fix);  

    
  } else
    throw std::runtime_error("Invalid input ID " + arrname +
			     std::string(" in compute/ft/spherical/bin."));


  Nstores = freq/every+1;


  // number of samples to be averaged over in each of the output values
  for (int i = 0; i < repeat; i++) 
    Nsamples.push_back(Nstores - i);

  
}


void FixCorrelateVector::setup()
{

  if (compute != nullptr) {

    input_array = compute->array.data();
    input_nc = compute->numberofcomponents;
    input_array_size = compute->array.size();
    

  } else if (fix != nullptr) {
    input_nc = fix->numberofcomponents;
    input_array_size = fix->array.size();
    input_array = fix->array.data();

  }
  

  // array stores the output for time auto-correlation 
  array.resize(repeat);

  for (auto &element : array)
    element = 0.0;


  // circular buffer for storage array
  storage_array.resize(repeat*input_array_size);
  circular_index = 0;  
  elements_in_storage = 0;

}


void FixCorrelateVector::start_of_step()
{
  
  if (integrate->timestep % every == 0) {
    if (compute != nullptr) {
      compute->this_step = true;
      compute->start_of_step();
    } else if (fix != nullptr) {
      fix->this_step = true;
      fix->start_of_step();
    }
  }

  if ((integrate->timestep - every) % freq == 0)  {

    for (auto &item : array)
      item = 0.0;

    // number of storage_array blocks which are storing useful data
    elements_in_storage = 1; 

    if (circular_index > 0) {
      double *point = get_storage_section(circular_index-1);
    
      double tmpsum = 0;
      for (int i = 0; i < input_array_size; i++)
	tmpsum += point[i]*point[i];

      array.at(0) = tmpsum;
    }
    
    
  }


  
  
}

double * FixCorrelateVector::get_storage_section(int index)
{
  return storage_array.data()+ (index % repeat)*input_array_size;
}


void FixCorrelateVector::add_to_array(int index,double *point1,
				      double *point2)
{
  double tmpsum = 0;
  
  for (int i = 0; i < input_array_size; i++)
    tmpsum += point1[i]*point2[i];
    
  array.at(index) += tmpsum;

}

void FixCorrelateVector::end_of_step()
{


  if (integrate->timestep % every != 0) return;


  double *point1 = get_storage_section(circular_index);
  // store the necessary info in the storage_array at the time specified
  for (int i = 0; i < input_array_size; i++)
    point1[i] = input_array[i];


  if (elements_in_storage < repeat)
    elements_in_storage += 1;
  
  
  for (int gamma = 0; gamma < elements_in_storage; gamma++) {

    double *point2 = get_storage_section(circular_index-gamma);

    add_to_array(gamma,point1,point2);
    
  }


  if (integrate->timestep % freq == 0) {
      for (int i = 0; i < repeat; i++)
        array.at(i) /= Nsamples.at(i);
  }
  
  if (compute != nullptr)	
    compute->this_step = false;
  else if (fix != nullptr)
    fix->this_step = false;
  
  
  circular_index += 1;  

}
