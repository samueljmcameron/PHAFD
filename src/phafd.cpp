
#include "phafd.hpp"
#include "atom.hpp"
#include "domain.hpp"
#include "grid.hpp"
#include "comm_brick.hpp"
#include "neighbor.hpp"
#include "input.hpp"
#include "integrate.hpp"
#include "group.hpp"
#include "fix.hpp"
#include "pair.hpp"
#include "compute.hpp"
#include "dump.hpp"

#include <algorithm>

using namespace PHAFD_NS;


PHAFD::PHAFD(MPI_Comm communicator, int argc, char **argv)
{
  world = communicator;

  inputfile = std::make_unique<std::fstream>();
  varmap = std::make_unique<std::map<std::string,std::string>>();



  int iarg = 1;  
  while(iarg < argc) {


    if (strcmp(argv[iarg],"-in") == 0) {
      if (iarg+1 == argc) {
	throw std::invalid_argument("Error: input flag '-in' specified, but no file given.");

      }
      inputfile->open(argv[iarg+1]);
      iarg += 2;
      
    } else if (strcmp(argv[iarg],"-var") == 0) {


      if (iarg + 2 >= argc) {
	throw std::invalid_argument("Error: invalid command line variable specification.");
      }


      std::string strVar = argv[iarg+2];


      int iiarg = 0;
      // check if double quotes (meaning variable with spaces in it)
      if (strVar.at(0) == '\"') {
	strVar.erase(strVar.begin());

	if (strVar.back() == '\"') {
	  strVar.pop_back();
	  if (strVar.empty())
	    throw std::invalid_argument("Empty double quotes in variable.");
	} else {


	  iiarg = 1;

	  while (iiarg + iarg + 2 < argc) {

	    strVar += std::string(" ") + std::string(argv[iarg+2+iiarg]);
	    if (strVar.back() == '\"') {
	      strVar.pop_back();
	      break;
	    }
	    
	    iiarg += 1;
	  
	  }

	  if (iarg + 2+ iiarg >= argc) {
	    throw std::invalid_argument("Error: open double quotes in command line variable specification.");
	  }

	  

	}
      }

      (*varmap)[argv[iarg+1]] = strVar;
      iarg += 3+iiarg;
    } else {
      throw std::invalid_argument("Error: invalid command line variable specification.");
    }
  }

  
  if (not inputfile->is_open()) {
    throw std::runtime_error("Error: need to specify input file.");
  }

  
  
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
  integrate = std::make_unique<Integrate>(this);

  return;
}

/*void PHAFD::init()
{

}
*/
