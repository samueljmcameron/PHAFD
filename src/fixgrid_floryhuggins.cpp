#include <algorithm>

#include "utility.hpp"


#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"

#include "fixgrid_floryhuggins.hpp"
#include "fftw_arr/array3d.hpp"

using namespace PHAFD_NS;


FixGridFloryHuggins::FixGridFloryHuggins(PHAFD *phafd) : Fix(phafd) {

};



void FixGridFloryHuggins::init(const std::vector<std::string> &v_line)
{


  Fix::init(v_line);

  std::vector<std::string> new_v_line(v_line);

  new_v_line.erase(new_v_line.begin());


  num_less_zero = num_great_zero = 0;
  temp = volFH = chi = 1;

  int iarg = 0;

  while (iarg < new_v_line.size()) {

    if (new_v_line[iarg] == "temp") {
      temp = std::stod(new_v_line[iarg+1]);
      iarg += 2;

    } else if (new_v_line[iarg] == "volFH") {
      volFH = std::stod(new_v_line[iarg+1]);
      iarg += 2;
    } else if (new_v_line[iarg] == "chi") {
      chi = std::stod(new_v_line[iarg+1]);
      iarg += 2;
    } else {
      throw std::runtime_error("Error: invalid fixgrid/floryhuggins command");
    }
  }


}

  
void FixGridFloryHuggins::setup()
{

}

  
void FixGridFloryHuggins::post_force()
{


  
  
  for (int i = 0; i < grid->chempot->Nz(); i++) 
    for (int j = 0; j < grid->chempot->Ny(); j++)
      for (int k = 0; k < grid->chempot->Nx(); k++) {
	if ((*grid->phi)(i,j,k) <= 0 ) {
	  (*grid->phi)(i,j,k) = -(*grid->phi)(i,j,k);
	  num_less_zero += 1;
	} else if ((*grid->phi)(i,j,k) >= 1) {
	  (*grid->phi)(i,j,k) = 1-((*grid->phi)(i,j,k)-1);
	  num_great_zero += 1;
	}
	(*grid->chempot)(i,j,k) 
	  += temp/volFH*(log((*grid->phi)(i,j,k)/(1-(*grid->phi)(i,j,k)))+chi*(1-2*(*grid->phi)(i,j,k)));
      }

}
