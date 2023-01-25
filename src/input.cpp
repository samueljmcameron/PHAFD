#include "input.hpp"


#include "utility.hpp"


#include "domain.hpp"
#include "grid.hpp"
#include "group.hpp"
#include "comm_brick.hpp"
#include "read_atoms.hpp"
#include "atom.hpp"
#include "neighbor.hpp"
#include "fixatom_semiflexible.hpp"
#include "fixgrid_conjugate.hpp"
#include "beadrodpmer/no_tether.hpp"
#include "ps_pde/conjugate_volfrac.hpp"
#include "pair_lj_cut.hpp"
#include "pair_gridatom_gaussian.hpp"
#include "pair_grid_floryhuggins.hpp"
#include "ps_pde/fixgrid_floryhuggins.hpp"

#include "ps_pde/iovtk.hpp"


#include <iostream>
#include <string>
#include <set>
#include <chrono>

using namespace PHAFD_NS;

Input::Input(PHAFD *phafd) : Pointers(phafd) {}

Input::~Input() = default;

void Input::read()
{
  std::vector<std::string> v_line;
  std::string line;

  line = "boxdims 100.0 100.0 100.0 boxorigin -50.0 -50.0 -50.0";
  v_line = utility::split_line(line);
  domain->set_box(v_line);


  line = "64 64 64 concentration";
  v_line = utility::split_line(line);

  grid->create(v_line);
  domain->set_subbox();
  // domain is now completely set.

  
  line = "constant concentration 0.2 0.2 seed 891409";
  v_line = utility::split_line(line);
  grid->populate(v_line);

  // grid is now completely set.


  std::string atom_data_fname = "read_p%.data";

  utility::replacePercentages(atom_data_fname,commbrick->me);
  

  ReadAtoms read_atoms(phafd);

  read_atoms.read_file(atom_data_fname);

  // atoms are now completely set.



  // need to figure out force cutoffs before neighbor and commbrick
  line = "4.0";
  pairs.push_back(std::make_unique<PairGridAtomGaussian>(phafd));
  v_line = utility::split_line(line);
  pairs.at(0)->settings(v_line);
  
  line = "0 phi 0.9 epsilon 4.0";
  v_line = utility::split_line(line);
  pairs.at(0)->coeff(v_line);

  line = "1 phi 0.1 epsilon 4.0";
  v_line = utility::split_line(line);
  pairs.at(0)->coeff(v_line);


  line = "2.0 1.0 1.781794";
  pairs.push_back(std::make_unique<PairLJCut>(phafd));
  v_line = utility::split_line(line);
  pairs.at(1)->settings(v_line);


  
  line = "temp 1.0 chi 3.0 volFH 0.01";
  pairs.push_back(std::make_unique<PairGridFloryHuggins>(phafd));
  v_line = utility::split_line(line);
  pairs.at(2)->settings(v_line);


  
  

  // set neighbor skin , infer cutoff from force cutoffs,
  //  and give pairs correct neigh_lists
  line = "2.0";
  v_line = utility::split_line(line);
  neighbor->setup(v_line);


  
  // can set commbrick now that neighbor cutoff is set
  commbrick->setup();



  // could set fixes and groups before or after neighbor and commbrick
  groups.push_back(std::make_unique<Group>(phafd));
  groups.at(0)->create_all();

  groups.push_back(std::make_unique<Group>(phafd));
  line = "nt_pmer0 molecule 0 1 2 3";
  v_line = utility::split_line(line);    
  groups.at(1)->create_group(v_line);


  groups.push_back(std::make_unique<Group>(phafd));
  line = "nt_pmer1 molecule 1";
  v_line = utility::split_line(line);    
  groups[2]->create_group(v_line);

  groups.push_back(std::make_unique<Group>(phafd));
  line = "nt_pmer2 molecule 2";
  v_line = utility::split_line(line);    
  groups[3]->create_group(v_line);

  groups.push_back(std::make_unique<Group>(phafd));
  line = "nt_pmer3 molecule 3";
  v_line = utility::split_line(line);    
  groups[4]->create_group(v_line);

  
  line = "fix_nt_pmer0 nt_pmer0 bondlength 2.0 zeta_para 0.5 zeta_perp 1.0 bending 206.0 temp 4.114 seed 489210";
  atomfixes.push_back(std::make_unique<FixAtomSemiFlexible<BeadRodPmer::NoTether>>(phafd));
  v_line = utility::split_line(line);
  atomfixes.at(0)->init(v_line);

  line = "fix_volfrac mobility 0.01 temp 1.0 volFH 0.01 gamma 10.0 seed 491240";
  gridfixes.push_back(std::make_unique<FixGridConjugate<psPDE::ConjugateVolFrac>>(phafd));
  v_line = utility::split_line(line);
  gridfixes.at(0)->init(v_line);
  

  
  // need to run things now?

  int step = 0;
  int Nsteps = 10000;
  double t = 0;
  double dt = 1e-3;

  
  domain->pbc();
  commbrick->borders();
  neighbor->build(step);


  
  // and save it
  std::string atom_prefix = "vtkfiles/atom_p%";
    
  utility::replacePercentages(atom_prefix,commbrick->me);

  std::string atom_fname = atom_prefix + std::string("_0.vtp");
  std::string atom_cname = atom_prefix + std::string(".pvd");

  atoms->Fs.setZero();
  psPDE::ioVTK::writeVTKcollectionHeader(atom_cname);
  writeVTKAtomData(atom_fname,"owned");
  psPDE::ioVTK::writeVTKcollectionMiddle(atom_cname,atom_fname,t);


  // and save it
  std::string field_prefix = "vtkfiles/field_p%";
  
  utility::replacePercentages(field_prefix,commbrick->me);
  
  std::string field_fname = field_prefix + std::string("_0.vti");
  std::string field_cname = field_prefix + std::string(".pvd");
  
  
  psPDE::ioVTK::writeVTKcollectionHeader(field_cname);
  psPDE::ioVTK::writeVTKImageData(field_fname,{grid->phi},domain->ps_domain->boxlo,
				  {grid->dx(),grid->dy(),grid->dz()});
  
  psPDE::ioVTK::writeVTKcollectionMiddle(field_cname,field_fname,0.0);
  


  for (auto &fix : atomfixes)
    fix->reset_dt(dt);

  for (auto &fix : atomfixes)
    fix->setup();


  for (auto &fix : gridfixes)
    fix->reset_dt(dt);

  for (auto &fix : gridfixes)
    fix->setup();

  if (commbrick->me == 0) 
    std::cout << "Running simulation of solution." << std::endl;


  int errflag = 0;
  int total_errflag;


  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  
  for (step = 1; step <= Nsteps; step++) {

    if (neighbor->decide()) {
      domain->pbc();
      commbrick->borders();
      neighbor->build(step);
    } else {
      commbrick->forward_comm();
    }


    atoms->Fs.setZero();

    for (auto & pair : pairs)
      pair->compute();

    commbrick->reverse_comm();

    for (auto &fix : atomfixes)
      fix->initial_integrate();

    grid->nonlinear->setZero();
    atoms->Fs.setZero();
    
    for (auto & pair : pairs)
      pair->compute();

    commbrick->reverse_comm();

    for (auto &fix : atomfixes)
      fix->final_integrate();

    for (auto &fix : gridfixes)
      fix->final_integrate();

    if (std::isnan((*grid->phi)(0,0,0)))
      errflag = 1;
	
    MPI_Allreduce(&errflag, &total_errflag,1,MPI_INT,MPI_SUM,world);

    if (total_errflag)
      throw std::runtime_error("NAN encountered in phi.");

    t += dt;





    if (step % 1000 == 0) {
      
      if (commbrick->me == 0)
	  std::cout << "Saving polymer on step " << step << std::endl;
      atoms->Fs.setZero();
      pairs.at(0)->compute();
      commbrick->reverse_comm();
      
      atom_fname = atom_prefix + std::string("_") + std::to_string(step) + std::string(".vtp");
      
      writeVTKAtomData(atom_fname,"owned");
      psPDE::ioVTK::writeVTKcollectionMiddle(atom_cname,atom_fname,t);
    }
    if (step % 1000 == 0) {
      if (commbrick->me == 0)
	std::cout << "Saving phi on step " << step << std::endl;
      
      field_fname = field_prefix + std::string("_") + std::to_string(step) + std::string(".vti");
      
      psPDE::ioVTK::writeVTKImageData(field_fname,{grid->phi},domain->ps_domain->boxlo,
				      {grid->dx(),grid->dy(),grid->dz()});
      
      psPDE::ioVTK::writeVTKcollectionMiddle(field_cname,field_fname,t);
      
    }	
    
    
    
  }

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

  
  psPDE::ioVTK::writeVTKcollectionFooter(atom_cname);    
  psPDE::ioVTK::writeVTKcollectionFooter(field_cname);

  if (commbrick->me == 0) {
    std::cout << "number of neighbor calls = " << neighbor->get_ncalls() << std::endl;
    
    std::cout << "Run time = "
	      << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1e6
	      << "seconds." << std::endl;  
  }


}


void Input::writeVTKAtomData(std::string fname,
			     std::string which_atoms)
/*============================================================================*/
/*
  Write scalar image data to a vtk (and paraview) compatible file
  (extension .vti).

  Parameters
  ----------

  fname : string
      Name of file to save with extension (either ".vtp" or ".pvtp").

  xs : beads to save
*/
/*============================================================================*/

{
  
  auto myfile = std::fstream(fname, std::ios::out);

  if (not myfile.is_open()) {
    throw std::runtime_error(std::string("Cannot open file ") + fname);
  }


  int numpoints;
  if (which_atoms == "all")
    numpoints = atoms->xs.cols();
  else if (which_atoms == "owned")
    numpoints = atoms->nowned;
  else if (which_atoms == "local")
    numpoints = atoms->nlocal;
  else if (which_atoms == "ghost")
    numpoints = atoms->nghost;
  else
    throw std::runtime_error("need which_atoms to be either all, owned, local, or ghost");
    
  
  myfile << "<?xml version=\"1.0\"?>" << std::endl
	 << "<VTKFile type=\"PolyData\"  version=\"1.0\""
	 << " byte_order=\"LittleEndian\">" << std::endl
	 << "<PolyData>" << std::endl
	 << "<Piece NumberOfPoints=\"" << numpoints
	 << "\" NumberOfLines=\"1\">" << std::endl;
  
  myfile << "<PointData Vectors=\"Force\">" << std::endl
	 << "<DataArray Name=\"Force\" "
	 << "type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
  
  if (which_atoms == "all") {
    for (int i = 0; i < atoms->Fs.cols(); i++) {
      myfile << atoms->Fs.col(i)(0) << " " << atoms->Fs.col(i)(1) << " "
	     << atoms->Fs.col(i)(2) << " ";
      
    }
  } else if (which_atoms == "owned") {
    for (int i = 0; i < atoms->nowned; i++) {
      if (atoms->labels[i] == PHAFD_NS::Atom::OWNED)
	myfile << atoms->Fs.col(i)(0) << " " << atoms->Fs.col(i)(1) << " "
	       << atoms->Fs.col(i)(2) << " ";
      else std::cout << "error?!" << std::endl;
    }
  } else if (which_atoms == "local") {
    for (int i = atoms->nowned; i < atoms->nowned+atoms->ngathered; i++) {
      if (atoms->labels[i] == PHAFD_NS::Atom::LOCAL)
	myfile << atoms->Fs.col(i)(0) << " " << atoms->Fs.col(i)(1) << " "
	       << atoms->Fs.col(i)(2) << " ";
    }
  } else if (which_atoms == "ghost") {
    for (int i = atoms->nowned; i < atoms->nowned+atoms->ngathered; i++) {
      if (atoms->labels[i] == PHAFD_NS::Atom::GHOST)
	myfile << atoms->Fs.col(i)(0) << " " << atoms->Fs.col(i)(1) << " "
	       << atoms->Fs.col(i)(2) << " ";
    } 
  }
  
  myfile << std::endl << "</DataArray>" << std::endl
	 << "</PointData>" << std::endl;

  myfile << "<Points>" << std::endl
	 << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"
	 << std::endl;

  if (which_atoms == "all") {
    for (int i = 0; i < atoms->xs.cols(); i++) {
      myfile << atoms->xs.col(i)(0) << " " << atoms->xs.col(i)(1) << " "
	     << atoms->xs.col(i)(2) << " ";
    
    }
  } else if (which_atoms == "owned") {
    for (int i = 0; i < atoms->nowned; i++) {
      if (atoms->labels[i] == PHAFD_NS::Atom::OWNED)
	myfile << atoms->xs.col(i)(0) << " " << atoms->xs.col(i)(1) << " "
	       << atoms->xs.col(i)(2) << " ";
      else std::cout << "error?!" << std::endl;
    }
  } else if (which_atoms == "local") {
    for (int i = atoms->nowned; i < atoms->nowned+atoms->ngathered; i++) {
      if (atoms->labels[i] == PHAFD_NS::Atom::LOCAL)
	myfile << atoms->xs.col(i)(0) << " " << atoms->xs.col(i)(1) << " "
	       << atoms->xs.col(i)(2) << " ";
    }
  } else if (which_atoms == "ghost") {
    for (int i = atoms->nowned; i < atoms->nowned+atoms->ngathered; i++) {
      if (atoms->labels[i] == PHAFD_NS::Atom::GHOST)
	myfile << atoms->xs.col(i)(0) << " " << atoms->xs.col(i)(1) << " "
	       << atoms->xs.col(i)(2) << " ";
    } 
  }

  myfile << std::endl << "</DataArray>" << std::endl
	 << "</Points>" << std::endl;

  myfile << "</Piece>" << std::endl
	 << "</PolyData>" << std::endl
	 << "</VTKFile>" << std::endl;    
  myfile.close();
  
}
