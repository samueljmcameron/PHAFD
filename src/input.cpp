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
#include "fixgrid_floryhuggins.hpp"
#include "beadrodpmer/no_tether.hpp"
#include "ps_pde/conjugate_volfrac.hpp"
#include "pair_lj_cut.hpp"
#include "pair_gridatom_gaussian.hpp"

#include "ps_pde/fixgrid_floryhuggins.hpp"

#include "ps_pde/iovtk.hpp"
#include "dump.hpp"

#include "compute_complex.hpp"
#include "fixgrid_ave.hpp"

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

  
  line = "constant concentration 0.2 0.2 891409";
  v_line = utility::split_line(line);
  grid->populate(v_line);

  // grid is now completely set.


  line = "atoms_p%.data";

  utility::replacePercentages(line,commbrick->me);
  

  ReadAtoms read_atoms(phafd);

  int errflag = read_atoms.read_file(line);
  if (errflag != ReadAtoms::SUCCESS) errflag = 1;
  else errflag = 0;
  int total_errflag;
  MPI_Allreduce(&errflag,&total_errflag,1,MPI_INT,MPI_SUM,world);
  if (total_errflag)
    throw std::runtime_error("Could not read atom file. ");



  atoms->Fs.setZero();
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
  line = "nt_pmers molecule 0 1 2";
  v_line = utility::split_line(line);    
  groups.at(1)->create_group(v_line);


  
  line = "fix_nt_pmers nt_pmers 489210 bondlength 2.0 zeta_para 0.5 zeta_perp 1.0 bending 206.0 temp 4.114";
  fixes.push_back(std::make_unique<FixAtomSemiFlexible<BeadRodPmer::NoTether>>(phafd));
  v_line = utility::split_line(line);
  fixes.at(0)->init(v_line);

  
  line = "fix_volfrac 491240 mobility 0.01 temp 1.0 volFH 0.01 gamma 10.0";
  fixes.push_back(std::make_unique<FixGridConjugate<psPDE::ConjugateVolFrac>>(phafd));
  v_line = utility::split_line(line);
  fixes.at(1)->init(v_line);


  

  line = "fix_chempot temp 1.0 chi 3.0 volFH 0.01";
  fixes.push_back(std::make_unique<FixGridFloryHuggins>(phafd));
  v_line = utility::split_line(line);
  fixes.at(2)->init(v_line);


  




  line = "atoms atom/owned vtkfiles/atom_p% 1000 attributes 4 x ux ix F";
  utility::replacePercentages(line,commbrick->me);

  dumps.push_back(std::make_unique<Dump>(phafd));
  v_line = utility::split_line(line);
  dumps.back()->init(v_line);



  line = "realgrid grid vtkfiles/field_p% 1000 attributes 2 phi nonlinear";
  utility::replacePercentages(line,commbrick->me);

  dumps.push_back(std::make_unique<Dump>(phafd));
  v_line = utility::split_line(line);
  dumps.back()->init(v_line);




  line = "fouriergrid ftgrid vtkfiles/ft_field_p% 1000 attributes 1 f_averagenorm";
  utility::replacePercentages(line,commbrick->me);

  dumps.push_back(std::make_unique<Dump>(phafd));
  v_line = utility::split_line(line);
  dumps.back()->init(v_line);

  
  line = "norm ft_phi norm";

  computes.push_back(std::make_unique<ComputeComplex>(phafd));
  v_line = utility::split_line(line);
  computes.back()->init(v_line);
  
  line = "averagenorm 100 10 1000 c_norm";
  fixes.push_back(std::make_unique<FixGridAve>(phafd));
  v_line = utility::split_line(line);
  fixes.back()->init(v_line);

  int step = 0;
  int Nsteps = 10000;
  double t = 0;
  double dt = 1e-3;
  

  // need to run things now?
  domain->pbc();
  commbrick->borders();
  neighbor->build(step);

  
  for (auto &dump : dumps) {
    dump->setup();
    dump->write_collection_header();
  }

  
  for (auto &fix : fixes)
    fix->reset_dt(dt);

  for (auto &fix : fixes)
    fix->setup();


  if (commbrick->me == 0) 
    std::cout << "Running simulation of solution." << std::endl;
  
    
  
  errflag = 0;
  total_errflag = 0;


  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();



  for (step = 1; step <= Nsteps; step++) {

    for (auto &compute: computes)
      compute->start_of_step();
    
    for (auto &fix: fixes)
      fix->start_of_step();


    for (auto &dump : dumps)
      dump->start_of_step(step);

    
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

    for (auto &fix : fixes)
      fix->initial_integrate();

    grid->nonlinear->setZero();
    atoms->Fs.setZero();
    
    for (auto & pair : pairs)
      pair->compute();

    commbrick->reverse_comm();

    // additional terms from e.g. chemical potential, which are added to grid->nonlinear
    for (auto &fix: fixes)
      fix->post_force();


    // mainly just fourier transforming phi and nonlinear
    for (auto &fix : fixes)
      fix->pre_final_integrate();


    // computes which act on fourier space grids
    for (auto &compute : computes)
      compute->in_fourier();

    
    // updating stuff (including fourier transformed phi and nonlinear
    for (auto &fix : fixes)
      fix->final_integrate();


    // mainly just inverse fourier transforming
    for (auto &fix : fixes)
      fix->post_final_integrate();

    for (auto &fix : fixes)
      fix->end_of_step();
    

    t += dt;

    
    for (auto &dump :dumps)
      if (step % dump->every == 0) {
	if (commbrick->me == 0)
	  std::cout << "Saving on step " << step << std::endl;


	dump->write_collection_middle(step);

      
      }
    
    if (std::isnan((*grid->phi)(0,0,0)))
      errflag = 1;
	
    MPI_Allreduce(&errflag, &total_errflag,1,MPI_INT,MPI_SUM,world);

    if (total_errflag) {

      for (auto &dump : dumps)
	dump->write_collection_footer();
      
      throw std::runtime_error("NAN encountered in phi.");
    }    
    
  }

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

  for (auto &dump :dumps)
    dump->write_collection_footer();
  
  if (commbrick->me == 0) {
    std::cout << "number of neighbor calls = " << neighbor->get_ncalls() << std::endl;
    
    std::cout << "Run time = "
	      << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1e6
	      << "seconds." << std::endl;  
  }


}
