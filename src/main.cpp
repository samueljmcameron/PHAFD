#include <vector>
#include <iostream>
#include <random>
#include <fstream>
#include <string>
#include <mpi.h>
#include <chrono>

#include "beadrodpmer/no_tether.hpp"
#include "beadrodpmer/single_tether.hpp"
#include "beadrodpmer/double_tether.hpp"

#include "ps_pde/grid.hpp"
#include "ps_pde/iovtk.hpp"
#include "ps_pde/conjugate_volfrac.hpp"
#include "ps_pde/fixgrid_floryhuggins.hpp"

#include "input.hpp"
#include "atom.hpp"
#include "comm_brick.hpp"
#include "domain.hpp"
#include "read_atoms.hpp"
#include "group.hpp"
#include "fixatom_semiflexible.hpp"
#include "neighbor.hpp"
#include "nbin.hpp"
#include "nstencil.hpp"
#include "neigh_list.hpp"
#include "npair.hpp"
#include "pair_lj_cut.hpp"
#include "pair_gridatom_gaussian.hpp"

void print_atoms(const Eigen::Matrix<double,3,Eigen::Dynamic> &,int);
void save_atoms(std::string ,const PHAFD::Atom &,int ,
		std::string which_atoms="all");

void writeVTKAtomData(std::string fname,
		      const PHAFD::Atom &atoms,
		      std::string which_atoms);

void make_vtk_name(std::string &fname ,const std::string & base,int step)
{
  fname = base + std::string("_") + std::to_string(step) + std::string(".vtp");
}


int main(int argc, char **argv)
{
  // imagine that the system is split into mpi_size rectangles

  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Info info = MPI_INFO_NULL;

  int ierr = MPI_Init(NULL,NULL);
  int mpi_size,id;
  ierr = MPI_Comm_size(comm,&mpi_size);
  ierr = MPI_Comm_rank(comm,&id);

  fftw_mpi_init();

  {
    // first set up the grid stuff
    
    std::string line = "boxdims 100.0 100.0 100.0 boxorigin -50.0 -50.0 -50.0";

    PHAFD::Domain domain(id,mpi_size,line);


    // then set up the atoms stuff

    std::unique_ptr<PHAFD::Atom> atoms;

    std::string atom_data_fname = "read_p%.data";

    PHAFD::input::replacePercentages(atom_data_fname,id);

    PHAFD::ReadAtoms read_atoms;


    read_atoms.read_file(atom_data_fname,atoms,comm,domain);

    atoms->final_init(comm,id,mpi_size);

    std::vector<PHAFD::Group> groups;

    groups.push_back(PHAFD::Group(*atoms));

    //    line = "points atom <> 110 149";

    //    groups.push_back(PHAFD::Group(line,*atoms));
	


    line = "nt_pmers molecule 0";

    groups.push_back(PHAFD::Group(line,*atoms));
    
    line = "nt_pmers molecule 1";

    groups.push_back(PHAFD::Group(line,*atoms));

    
    line = "nt_pmers molecule 2";

    groups.push_back(PHAFD::Group(line,*atoms));


    line = "nt_pmers molecule 3";

    groups.push_back(PHAFD::Group(line,*atoms));



    std::vector<std::unique_ptr<PHAFD::FixAtom>> fixatoms;


    int seed0 = 124980;
    int seed1 = seed0 + 3910;
    int seed2 = seed1 + 3910;
    int seed3 = seed2 + 3910;

    
    line = "bondlength 2.0 zeta_para 0.5 zeta_perp 1.0 bending 206.0 temp 4.114 seed "
      + std::to_string(seed0);
    
    fixatoms.push_back(std::make_unique<PHAFD::FixAtomSemiFlexible
		       <BeadRodPmer::NoTether>>(groups[1],line,*atoms,domain));

    line = "bondlength 2.0 zeta_para 0.5 zeta_perp 1.0 bending 206.0 temp 4.114 seed "
      + std::to_string(seed1);

    fixatoms.push_back(std::make_unique<PHAFD::FixAtomSemiFlexible
		       <BeadRodPmer::NoTether>>(groups[2],line,*atoms,domain));


    line = "bondlength 2.0 zeta_para 0.5 zeta_perp 1.0 bending 206.0 temp 4.114 seed "
      + std::to_string(seed2);

    fixatoms.push_back(std::make_unique<PHAFD::FixAtomSemiFlexible
		       <BeadRodPmer::NoTether>>(groups[3],line,*atoms,domain));


    line = "bondlength 2.0 zeta_para 0.5 zeta_perp 1.0 bending 206.0 temp 4.114 seed "
      + std::to_string(seed3);

    fixatoms.push_back(std::make_unique<PHAFD::FixAtomSemiFlexible
		       <BeadRodPmer::NoTether>>(groups[4],line,*atoms,domain));


    // apply pbc right away
    domain.pbc(*atoms);

    // and save it
    std::string atom_prefix = "vtkfiles_continuation/atom_p%";
    
    PHAFD::input::replacePercentages(atom_prefix,id);

    std::string atom_fname = atom_prefix + std::string("_0.vtp");
    std::string atom_cname = atom_prefix + std::string(".pvd");


    atoms->Fs.setZero();
    psPDE::ioVTK::writeVTKcollectionHeader(atom_cname);
    writeVTKAtomData(atom_fname,*atoms,"owned");
    psPDE::ioVTK::writeVTKcollectionMiddle(atom_cname,atom_fname,0.0);



    // then set up grid stuff

    line = "64 64 64 concentration";

    std::vector<std::string> v_line = PHAFD::input::split_line(line);

    psPDE::Grid grid(v_line,comm);

    
    
    line = "read concentration vtkfiles/field_p%_200000.vti";

    PHAFD::input::replacePercentages(line,id);

    
    v_line = PHAFD::input::split_line(line);

    
    
    grid.populate(v_line,domain);
    domain.partition(grid);

    // and save it
    std::string field_prefix = "vtkfiles_continuation/field_p%";
    
    PHAFD::input::replacePercentages(field_prefix,id);

    std::string field_fname = field_prefix + std::string("_0.vti");
    std::string field_cname = field_prefix + std::string(".pvd");


    psPDE::ioVTK::writeVTKcollectionHeader(field_cname);
    psPDE::ioVTK::writeVTKImageData(field_fname,{grid.phi.get()},domain.boxlo,
				    {grid.dx(),grid.dy(),grid.dz()});

    psPDE::ioVTK::writeVTKcollectionMiddle(field_cname,field_fname,0.0);



    psPDE::ConjugateVolFrac conjvfrac(domain,grid);


    line = "mobility 0.01 temp 1.0 volFH 0.01 gamma 10.0 seed 491240";

    v_line = PHAFD::input::split_line(line);

    conjvfrac.readCoeffs(v_line);





    psPDE::FixGridFloryHuggins fxgridFH;

    line = "temp 1.0 chi 3.0 volFH 0.01";

    v_line = PHAFD::input::split_line(line);

    fxgridFH.readCoeffs(v_line);

    double dt = 1e-3;
    
    double t = 0;


    // define interaction cut off radius
    
    double globalcut = 4.0;
    
    double ljcut = 2.0;
    double skin = 2.0;
    double sigma = ljcut/1.1226;
    
    
    double gacut = globalcut;
    
    domain.pbc(*atoms);
    
    
    // set up inter-processor communication details
    PHAFD::CommBrick commbrick(comm);
    commbrick.setup(domain,globalcut);
    commbrick.borders(*atoms);
    
    
    PHAFD::Neighbor neighbor;
    neighbor.setup(domain,globalcut,skin);
    neighbor.build(*atoms,grid,0);  
    
    
    
    // read in lj parameters and grid atom parameters
    std::vector<std::string> pairljsettings =  {"2.0","1.0","1.781794"};
    std::vector<std::string> pairgridatomsettings = {"4.0"};
    
    
    std::vector<std::string> pairgridatomcoeff0 = {"0","phi","0.9","epsilon","4.0"};
    std::vector<std::string> pairgridatomcoeff1 = {"1","phi","0.1","epsilon","4.0"};
    std::vector<std::string> pairgridatomcoeff2 = {"2","epsilon","0.0"};
    
    
    
    
    // read in pair data
    PHAFD::PairGridAtomGaussian pairgridatom(atoms.get(),&grid);
    pairgridatom.settings(pairgridatomsettings);
    pairgridatom.coeff(pairgridatomcoeff0);
    pairgridatom.coeff(pairgridatomcoeff1);
    //    pairgridatom.coeff(pairgridatomcoeff2);
    pairgridatom.init_list(neighbor.neigh_lists[0].get());
    
    
    
    
    PHAFD::PairLJCut pairlj(atoms.get(),nullptr);
    pairlj.settings(pairljsettings);
    pairlj.init_list(neighbor.neigh_lists[1].get());
    
    
    
    int step = 0;
    int Nsteps = 200000;
    
    conjvfrac.reset_dt(dt);

    for (auto &fix : fixatoms)
      fix->reset_dt(dt);
    

    for (auto &fix : fixatoms)
      fix->setup();
  

    
    
    int errflag = 0;
    
    int totalerr;


    if (id == 0){
      std::cout << "Running simulation of solution." << std::endl;
    }
    
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    
    for (int step = 1; step <= Nsteps; step++) {


      if (neighbor.decide(*atoms,comm)) {
	domain.pbc(*atoms);
	commbrick.borders(*atoms);
	neighbor.build(*atoms,grid,step);
	
	
      } else {
	commbrick.forward_comm(*atoms);
      }

      // compute forces on atoms from configuration (both grid and atom) at time t
      atoms->Fs.setZero();
      pairlj.compute(domain);      
      pairgridatom.compute(domain);      
      commbrick.reverse_comm(*atoms);


      // then integrate to get first half step, R(t + 0.5*dt)
      for (auto &fix : fixatoms)
	fix->initial_integrate();


      grid.nonlinear->setZero();
      // compute forces on atoms from configuration (both grid and atom) at time t+0.5*dt
      atoms->Fs.setZero();
      pairlj.compute(domain);
      pairgridatom.compute(domain);
      commbrick.reverse_comm(*atoms);


      // then integrate to get second half step, R(t + dt)
      for (auto &fix : fixatoms)
	fix->final_integrate();      
      
      // add flory huggins contribution to nonlinear, which already has contribution at t+0.5*dt
      fxgridFH.compute(grid);      

      // FFT phi(r,t) 
      
      fftw_execute(grid.forward_phi);
      fftw_execute(grid.forward_nonlinear);


      // compute phi(q,t+dt) 

      conjvfrac.update();

      // IFFT back to get phi(r,t+dt) 
      fftw_execute(grid.backward_phi);
      
      if (std::isnan((*grid.phi)(0,0,0) )) 
	errflag = 1;

      MPI_Allreduce(&errflag,&totalerr,1,MPI_INT,MPI_SUM,grid.comm);

      if (totalerr) {
	psPDE::ioVTK::writeVTKcollectionFooter(field_cname);
	throw std::runtime_error("NAN encountered in phi.");

      }
    
      t += dt;


      if (step % 1000 == 0) {

	if (id == 0)
	  std::cout << "Saving polymer on step " << step << std::endl;
	atoms->Fs.setZero();
	pairgridatom.compute(domain);
	commbrick.reverse_comm(*atoms);
	
	atom_fname = atom_prefix + std::string("_") + std::to_string(step) + std::string(".vtp");

	writeVTKAtomData(atom_fname,*atoms,"owned");
	psPDE::ioVTK::writeVTKcollectionMiddle(atom_cname,atom_fname,t);
      }
      if (step % 1000 == 0) {
	if (id == 0)
	  std::cout << "Saving phi on step " << step << std::endl;

	field_fname = field_prefix + std::string("_") + std::to_string(step) + std::string(".vti");

	psPDE::ioVTK::writeVTKImageData(field_fname,{grid.phi.get()},domain.boxlo,
					{grid.dx(),grid.dy(),grid.dz()});
	
	psPDE::ioVTK::writeVTKcollectionMiddle(field_cname,field_fname,t);
	
      }	

    }
    psPDE::ioVTK::writeVTKcollectionFooter(atom_cname);    
    psPDE::ioVTK::writeVTKcollectionFooter(field_cname);

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    
    if (id == 0) {
      std::cout << "number of neighbor calls = " << neighbor.get_ncalls() << std::endl;

      std::cout << "Run time = "
		<< std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1e6
		<< "seconds." << std::endl;  
    }


    
    
  }
  fftw_mpi_cleanup();

  
  ierr = MPI_Finalize();
  return 0;

}



void writeVTKAtomData(std::string fname,
		      const PHAFD::Atom &atoms,
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
    numpoints = atoms.xs.cols();
  else if (which_atoms == "owned")
    numpoints = atoms.nowned;
  else if (which_atoms == "local")
    numpoints = atoms.nlocal;
  else if (which_atoms == "ghost")
    numpoints = atoms.nghost;
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
    for (int i = 0; i < atoms.Fs.cols(); i++) {
      myfile << atoms.Fs.col(i)(0) << " " << atoms.Fs.col(i)(1) << " "
	     << atoms.Fs.col(i)(2) << " ";
      
    }
  } else if (which_atoms == "owned") {
    for (int i = 0; i < atoms.nowned; i++) {
      if (atoms.labels[i] == PHAFD::Atom::OWNED)
	myfile << atoms.Fs.col(i)(0) << " " << atoms.Fs.col(i)(1) << " "
	       << atoms.Fs.col(i)(2) << " ";
      else std::cout << "error?!" << std::endl;
    }
  } else if (which_atoms == "local") {
    for (int i = atoms.nowned; i < atoms.nowned+atoms.ngathered; i++) {
      if (atoms.labels[i] == PHAFD::Atom::LOCAL)
	myfile << atoms.Fs.col(i)(0) << " " << atoms.Fs.col(i)(1) << " "
	       << atoms.Fs.col(i)(2) << " ";
    }
  } else if (which_atoms == "ghost") {
    for (int i = atoms.nowned; i < atoms.nowned+atoms.ngathered; i++) {
      if (atoms.labels[i] == PHAFD::Atom::GHOST)
	myfile << atoms.Fs.col(i)(0) << " " << atoms.Fs.col(i)(1) << " "
	       << atoms.Fs.col(i)(2) << " ";
    } 
  }
  
  myfile << std::endl << "</DataArray>" << std::endl
	 << "</PointData>" << std::endl;

  myfile << "<Points>" << std::endl
	 << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"
	 << std::endl;

  if (which_atoms == "all") {
    for (int i = 0; i < atoms.xs.cols(); i++) {
      myfile << atoms.xs.col(i)(0) << " " << atoms.xs.col(i)(1) << " "
	     << atoms.xs.col(i)(2) << " ";
    
    }
  } else if (which_atoms == "owned") {
    for (int i = 0; i < atoms.nowned; i++) {
      if (atoms.labels[i] == PHAFD::Atom::OWNED)
	myfile << atoms.xs.col(i)(0) << " " << atoms.xs.col(i)(1) << " "
	       << atoms.xs.col(i)(2) << " ";
      else std::cout << "error?!" << std::endl;
    }
  } else if (which_atoms == "local") {
    for (int i = atoms.nowned; i < atoms.nowned+atoms.ngathered; i++) {
      if (atoms.labels[i] == PHAFD::Atom::LOCAL)
	myfile << atoms.xs.col(i)(0) << " " << atoms.xs.col(i)(1) << " "
	       << atoms.xs.col(i)(2) << " ";
    }
  } else if (which_atoms == "ghost") {
    for (int i = atoms.nowned; i < atoms.nowned+atoms.ngathered; i++) {
      if (atoms.labels[i] == PHAFD::Atom::GHOST)
	myfile << atoms.xs.col(i)(0) << " " << atoms.xs.col(i)(1) << " "
	       << atoms.xs.col(i)(2) << " ";
    } 
  }

  myfile << std::endl << "</DataArray>" << std::endl
	 << "</Points>" << std::endl;

  myfile << "</Piece>" << std::endl
	 << "</PolyData>" << std::endl
	 << "</VTKFile>" << std::endl;    
  myfile.close();
  
}
