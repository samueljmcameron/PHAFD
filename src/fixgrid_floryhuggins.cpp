#include <algorithm>

#include "utility.hpp"


#include "grid.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"

#include "fixgrid_floryhuggins.hpp"
#include "ps_pde/fixgrid_floryhuggins.hpp"

using namespace PHAFD_NS;


FixGridFloryHuggins::FixGridFloryHuggins(PHAFD *phafd) : FixGrid(phafd) {

  ps_flory = std::make_unique<psPDE::FixGridFloryHuggins>();
};



void FixGridFloryHuggins::init(const std::vector<std::string> &v_line)
{


  FixGrid::init(v_line);

  std::vector<std::string> new_v_line(v_line);

  new_v_line.erase(new_v_line.begin());
  
  ps_flory->readCoeffs(new_v_line);

}

  
void FixGridFloryHuggins::setup()
{

}

  
void FixGridFloryHuggins::pre_final_integrate()
{

  ps_flory->compute(*(grid->ps_grid));
}
