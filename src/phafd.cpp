
#include "phafd.hpp"
#include "atom.hpp"
#include "domain.hpp"
#include "grid.hpp"
#include "comm_brick.hpp"
#include "neighbor.hpp"
#include "input.hpp"


using namespace PHAFD_NS;


PHAFD::PHAFD(MPI_Comm communicator)
{
  world = communicator;
  create();

}

PHAFD::~PHAFD()  = default;

void PHAFD::create()
{


  atoms = std::make_unique<Atom>(this);
  grid = std::make_unique<Grid>(this);
  commbrick = std::make_unique<CommBrick>(this);
  domain = std::make_unique<Domain>(this);
  neighbor = std::make_unique<Neighbor>(this);
  input = std::make_unique<Input>(this);

  return;
}

/*void PHAFD::init()
{

}
*/
