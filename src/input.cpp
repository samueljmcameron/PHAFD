#include "input.hpp"


#include "utility.hpp"


#include "domain.hpp"
#include "grid.hpp"
#include "group.hpp"
#include "comm_brick.hpp"
#include "atom.hpp"
#include "neighbor.hpp"
#include "fixatom_semiflexible.hpp"
#include "fixgrid_conjugate.hpp"
#include "fixgrid_floryhuggins.hpp"
#include "fixgrid_gradphi.hpp"
#include "fixatom_drag.hpp"
#include "beadrodpmer/no_tether.hpp"
#include "beadrodpmer/single_tether.hpp"
#include "beadrodpmer/double_tether.hpp"


#include "conjugate_volfrac.hpp"
#include "pair_lj_cut.hpp"
#include "pair_harmonic_cut.hpp"
#include "pair_gridatom_gaussian.hpp"
#include "pair_gridatom_lj_ish.hpp"
#include "pair_gridatom_lj_ish_linear.hpp"
#include "fixgrid_floryhuggins.hpp"

#include "dump.hpp"
#include "integrate.hpp"
#include "compute_complex.hpp"
#include "compute_pair.hpp"
#include "fixgrid_ave.hpp"

#include <string>
#include <set>
#include <chrono>

using namespace PHAFD_NS;

Input::Input(PHAFD *phafd) : Pointers(phafd) {}

Input::~Input() = default;

void Input::read()
{


  std::vector<std::string> v_line;
  std::string line, firstword;
  
  while(std::getline(*inputfile,line)) {
    
    utility::convertVariables(line,*varmap);
    utility::replacePercentages(line,commbrick->me);

    v_line = utility::split_line(line);

    if (v_line.size() == 0) continue;
    
    firstword = v_line.at(0);
    v_line.erase(v_line.begin());
    
    if (firstword == "domain") {
      domain->set_box(v_line);
    } else if (firstword == "grid_style") {
      if (!domain->boxset)
	throw std::invalid_argument("CAnnot use command grid_style before domain box is set.");
      grid->create(v_line);
      domain->set_subbox();
      
    } else if (firstword == "grid_populate") {
      if (!grid->gridset)
	throw std::invalid_argument("CAnnot use command grid_populate before grid_style is set.");


      grid->populate(v_line);
    } else if (firstword == "atom_style") {

      // set atom style and number of atoms on the processor.
      atoms->setup(v_line);
    } else if (firstword == "atom_populate") {

      if  (!atoms->atomset)
	throw std::invalid_argument("CAnnot use command atom_populate before atom_style is set.");


      atoms->populate(v_line);

      // call this to ensure everything is filled in that needs to be filled in.
      atoms->check_tags_and_types();

      if (commbrick->me == 0 && atoms->ntypes == 0)
	std::cout << "WARNING: running simulation with no atoms." << std::endl;
      
    } else if (firstword == "pair") {
      
      
      if (atoms->ntypes <= 0)
	throw std::runtime_error("Cannot create pair style before atoms are made.");
      
      std::string pairname = v_line.at(0);
      v_line.erase(v_line.begin());
      if (pairname == "gridatom/gaussian") {
	  
	pairs.push_back(std::make_unique<PairGridAtomGaussian>(phafd));
	
      } else if (pairname == "gridatom/LJish") {
	  
	pairs.push_back(std::make_unique<PairGridAtomLJish>(phafd));
	
      } else if (pairname == "gridatom/LJish/linear") {
	  
	pairs.push_back(std::make_unique<PairGridAtomLJishLinear>(phafd));
	
      } else if (pairname == "lj/cut") {
	
	pairs.push_back(std::make_unique<PairLJCut>(phafd));
	
      } else if (pairname == "harmonic/cut") {
	
	pairs.push_back(std::make_unique<PairHarmonicCut>(phafd));
	
      } else
	throw std::invalid_argument("Unrecognised pair style.");
      
      
      pairs.back()->settings(v_line);
      
      std::streampos fpos =  inputfile->tellg();
      while (std::getline(*inputfile,line)) {
	utility::convertVariables(line,*varmap);
	
	if (line == "" || line == "#") continue;
	v_line = utility::split_line(line);
	
	firstword = v_line.at(0);
	v_line.erase(v_line.begin());



	
	if (firstword == "pairCoeffs") {
	  pairs.back()->coeff(v_line);
	  fpos = inputfile->tellg();
	} else {
	  inputfile->seekg(fpos);
	  break;
	}
      }
    } else if (firstword == "neighbor") {
      if (pairs.size() == 0 && atoms->ntypes > 0)
	throw std::runtime_error("Must call neighbor after pairs are given.");
      
      if (!domain->subboxset)
	throw std::runtime_error("Must call neighbor after subboxes are set (via grid style).");
      
      neighbor->setup(v_line);
      // can set commbrick now that neighbor cutoff is set
      commbrick->setup();
      
    } else if (firstword == "group") {
      
      if (atoms->ntypes <= 0)
	throw std::runtime_error("Must call group after atoms have been read in (via read_atoms).");
      
      groups.push_back(std::make_unique<Group>(phafd));
      groups.back()->create_group(v_line);
      
	
      
    } else if (firstword == "fix") {
      
      firstword = v_line.at(0);
      v_line.erase(v_line.begin());
      
      if (firstword == "atom/semiflexible/notether") {
	if (atoms->ntypes <= 0)
	  throw std::runtime_error("Fix requires atoms to be created.");
	
	
	fixes.push_back(std::make_unique<FixAtomSemiFlexible<BeadRodPmer::NoTether>>(phafd));
	
      } else if (firstword == "atom/semiflexible/doubletether") {
	if (atoms->ntypes <= 0)
	  throw std::runtime_error("Fix requires atoms to be created.");
	
	
	fixes.push_back(std::make_unique<FixAtomSemiFlexible<BeadRodPmer::DoubleTether>>(phafd));
	
      } else if (firstword == "atom/semiflexible/singletether") {
	if (atoms->ntypes <= 0)
	  throw std::runtime_error("Fix requires atoms to be created.");
	
	
	fixes.push_back(std::make_unique<FixAtomSemiFlexible<BeadRodPmer::SingleTether>>(phafd));
	
      } else if (firstword == "atom/drag") {
	if (atoms->ntypes <= 0)
	  throw std::runtime_error("Fix requires atoms to be created.");
	
	
	fixes.push_back(std::make_unique<FixAtomDrag>(phafd));
	
      

      } else if (firstword == "grid/conjugate/volfrac") {
	if (!grid->gridpopulated)
	  throw std::runtime_error("Fix requires grid to be populated.");
	fixes.push_back(std::make_unique<FixGridConjugate<ConjugateVolFrac>>(phafd));
	
      } else if (firstword == "grid/floryhuggins") {
	if (!grid->gridpopulated)
	  throw std::runtime_error("Fix requires grid to be populated.");
	fixes.push_back(std::make_unique<FixGridFloryHuggins>(phafd));	  
      } else if (firstword == "grid/ave") {
	fixes.push_back(std::make_unique<FixGridAve>(phafd));
      } else if (firstword == "grid/gradphi") {
	if (!grid->gridpopulated)
	  throw std::runtime_error("Fix requires grid to be populated.");
	fixes.push_back(std::make_unique<FixGridGradPhi>(phafd));
      }	else
	throw std::runtime_error("Invalid fix.");
      
      
      fixes.back()->init(v_line);
      
    } else if (firstword == "dump") {
      dumps.push_back(std::make_unique<Dump>(phafd));
      dumps.back()->init(v_line);
    } else if (firstword == "compute") {
      
      firstword = v_line.at(0);
      v_line.erase(v_line.begin());
      
      if (firstword == "complex")
	
	computes.push_back(std::make_unique<ComputeComplex>(phafd));
      
      else if (firstword == "pair")
	
	computes.push_back(std::make_unique<ComputePair>(phafd));

      else
	throw std::runtime_error("Invalid compute.");
      
      
      computes.back()->init(v_line);
      
    } else if (firstword == "timestep") {
      
      integrate->firststep = std::stoi(v_line.at(0));
      integrate->timestep = integrate->firststep;
      
    } else if (firstword == "dt") {
      
      integrate->dt = std::stod(v_line.at(0));
      
    } else if (firstword == "run") {
      
      if (integrate->dt <= 0.0)
	throw std::runtime_error("Cannot run simulation without setting dt > 0.");
      
      integrate->nsteps = std::stoi(v_line.at(0));
      
      integrate->setup();
      integrate->run();
    } else if (firstword == "run_until_touching") {
      
      if (integrate->dt <= 0.0)
	throw std::runtime_error("Cannot run simulation without setting dt > 0.");
      
      integrate->nsteps = std::stoi(v_line.at(0));

      integrate->setup();
      integrate->run_until_touching(std::stod(v_line.at(1)));
    
    }
  }
  return;
}
