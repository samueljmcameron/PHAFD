#include <algorithm>

#include "utility.hpp"


#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "fix.hpp"
#include "compute.hpp"
#include "dump.hpp"

#include "fixgrid_ave.hpp"

using namespace PHAFD_NS;

FixGridAve::FixGridAve(PHAFD *phafd) : Fix(phafd) {
  averaging = true;
};


void FixGridAve::init(const std::vector<std::string> &v_line)
{

  compute = nullptr;
  fix = nullptr;

  Fix::init(v_line);

  std::vector<std::string> new_v_line(v_line);

  every = std::stoi(new_v_line.at(1));
  repeat = std::stoi(new_v_line.at(2));
  freq = std::stoi(new_v_line.at(3));

  if (every <= 0 || repeat <= 0 || freq <= 0)
    throw std::invalid_argument("Every, repeat, and freq must all be >= 0.");  
  if (every*repeat > freq)
    throw std::invalid_argument("Every and repeat exceed freq in fixgridave.");

  std::string quantity_name = new_v_line.at(4);

  if (quantity_name.rfind("c_",0)==0) {
      
    for (auto &other : computes)
      if (quantity_name.substr(2) == other->name)
	compute = other.get();

    if (compute == nullptr)
      throw std::invalid_argument("Cannot compute average of " + quantity_name +
				  std::string("as no matching compute ID exists."));

    if (compute->per_grid)
      per_grid = true;
    else if (compute->per_ftgrid)
      per_ftgrid = true;
    else
      throw std::invalid_argument("Need per_grid or per_ftgrid compute for fixgridave.");


    Nx = compute->Nx;
    Ny = compute->Ny;
    Nz = compute->Nz;


    
  } else if (quantity_name.rfind("f_",0) == 0) {
    for (auto &other : fixes)
      if (quantity_name.substr(2) == other->name)
	fix = other.get();

    if (fix == nullptr)
      throw std::invalid_argument("Cannot compute average of " + quantity_name +
				  std::string("as no matching fix ID exists."));

    
    if (fix->per_grid)
      per_grid = true;
    else if (fix->per_ftgrid)
      per_ftgrid = true;
    else
      throw std::invalid_argument("Need per_grid or per_ftgrid fix for fixgridave.");


    Nx = fix->Nx;
    Ny = fix->Ny;
    Nz = fix->Nz;
    
  }

  array.resize(Nx*Ny*Nz);




}

  
void FixGridAve::setup()
{



  for (auto &dump : dumps) {
    if (dump_callers.find(dump->name) != dump_callers.end())
      if (dump->every != freq)
	throw std::runtime_error("Using fixgrid/ave with frequency that doesn't match "
				 "dump frequency - waste of time.");
  }
  
  step_counter = 0;
  next_freq = freq;


  savesteps.clear();

  for (int i = repeat-1; i >= 0; i--)
    savesteps.push_back(next_freq-every*i);
  
  for (int i  = 0 ; i < array.size(); i++)
    array[i] = 0.0;


  collect_this_step = false;
  
}


void FixGridAve::start_of_step()
{

  step_counter += 1;
  
  collect_this_step = false;
  
  for (auto save : savesteps) {
    
    if (step_counter == save) {
      
      if (compute != nullptr)	
	compute->this_step = true;
      else if (fix != nullptr)
	fix->this_step = true;

      collect_this_step = true;

      break;
    }
  }


  
  
}


void FixGridAve::end_of_step() {




  if (step_counter == savesteps.at(0)) {

    for (int i  = 0 ; i < array.size(); i++)
      array[i] = 0.0;
    

  }

  if (step_counter == savesteps.back()) {
    next_freq += freq;

    savesteps.clear();
    for (int i = repeat-1; i >= 0; i--)
      savesteps.push_back(next_freq-every*i);


  
  
    if (compute != nullptr && collect_this_step)	 {
      
      for (int i  = 0 ; i < array.size(); i++) {
	array[i] += compute->array[i];
	array[i] /= repeat;
      }
      
    } else if (fix != nullptr && collect_this_step) {
      
      for (int i  = 0 ; i < array.size(); i++) {
	array[i] += fix->array[i];
	array[i] /= repeat;
      }
    }
    
  } else {

    if (compute != nullptr && collect_this_step)	 {
      
      for (int i  = 0 ; i < array.size(); i++) {
	array[i] += compute->array[i];

      }
      
    } else if (fix != nullptr && collect_this_step) {
      
      for (int i  = 0 ; i < array.size(); i++) {
	array[i] += fix->array[i];

      }
    }
  }

  

}
