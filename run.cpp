#include <fstream>
#include <memory>

#include "ps_pde/fftw_mpi_3darray.hpp"
#include "ps_pde/integrator.hpp"

#include "ps_pde/conjplane.hpp"
#include "ps_pde/randompll.hpp"
#include "ps_pde/griddata.hpp"
#include "ps_pde/iovtk.hpp"
#include "ps_pde/timestep.hpp"

#include "beadrodpmer/iovtk.hpp"


#include "poly_wrappers.hpp"

#include "run.hpp"


void transfer_nucleation_site(std::vector<double> & X_i,
			 const Eigen::Vector3d & pR,
			 double Lx, double Ly, double Lz)
{
  X_i[0] = pR(0);
  X_i[1] = pR(1);
  X_i[2] = pR(2);

  while (X_i[0] > Lx/2)
    X_i[0] -= Lx;
  while (X_i[0] < -Lx/2)
    X_i[0] += Lx;
  while (X_i[1] > Ly/2)
    X_i[1] -= Ly;
  while (X_i[1] < -Ly/2)
    X_i[1] += Ly;
  while (X_i[2] > Lz/2)
    X_i[2] -= Lz;
  while (X_i[2] < -Lz/2)
    X_i[2] += Lz;
}


void transfer_free_energy(std::vector<double> & dFdX,
			  const std::vector<std::vector<double>> & free_energy_derivs,
			  int index)
{

  dFdX[0] = free_energy_derivs[index][0];
  dFdX[1] = free_energy_derivs[index][1];
  dFdX[2] = free_energy_derivs[index][2];

}


void single_step(double t, NoTetherWrap & pmer,
		 double dt)
{
  pmer.set_unprojected_noise(dt);
  pmer.update_G();
  pmer.update_Hhat();
  pmer.compute_noise();
  pmer.compute_effective_kappa();
  pmer.compute_uc_forces();
  pmer.add_external_force(pmer.dFdX,pmer.nuc_beads[0]);
  
  pmer.compute_tension();
  pmer.initial_integrate(dt);
  

  
  pmer.update_Hhat();
  pmer.compute_effective_kappa();
  pmer.compute_uc_forces();
  pmer.add_external_force(pmer.dFdX,pmer.nuc_beads[0]);
  
  pmer.compute_tension();

  int itermax = 20;
  
  int iterations = pmer.correct_tension(dt,itermax,1e-8);

  
  if (iterations > itermax) {
    std::cout << "too many iterations when correcting tension at time " << t
	      <<  ", retrying the step with new noise " << std::endl;
    for (int i = 0; i < pmer.get_Nbeads(); i++) 
      pmer.atoms[i].R = pmer.Rtmp[i];

    single_step(t,pmer,dt);
  }
  else {
    pmer.final_integrate(dt);
  
    pmer.compute_tangents_and_friction();
  }

  return;
}

void single_step(double  t,DoubleTetherWrap & pmer,
		 double dt)
{
  pmer.set_unprojected_noise(dt);
  pmer.update_G();
  pmer.update_Hhat();
  pmer.compute_noise();
  pmer.compute_effective_kappa();
  pmer.compute_uc_forces();
  pmer.add_external_force(pmer.dFdX,pmer.nuc_beads[0]);
  
  
  pmer.compute_tension({0.0,0.0,0.0},{0.0,0.0,0.0});
  pmer.initial_integrate(dt);
  

  
  pmer.update_Hhat();
  pmer.compute_effective_kappa();
  pmer.compute_uc_forces();
  pmer.add_external_force(pmer.dFdX,pmer.nuc_beads[0]);
  
  pmer.compute_tension({0.0,0.0,0.0},{0.0,0.0,0.0});
  int itermax = 20;
  
  int iterations = pmer.correct_tension(dt,pmer.atoms[0].R,
					pmer.atoms[pmer.get_Nbeads()-1].R,
					itermax,1e-8);

  
  if (iterations > itermax) {
    std::cout << "too many iterations when correcting tension at time " << t
	      <<  ", retrying the step with new noise " << std::endl;
    for (int i = 0; i < pmer.get_Nbeads(); i++) 
      pmer.atoms[i].R = pmer.Rtmp[i];

    single_step(t,pmer,dt);
  }
  else {
    pmer.final_integrate(dt);
  
    pmer.compute_tangents_and_friction();
  }

  return;
}


void single_step(double t,SingleTetherWrap & pmer,
		 double dt)
{
  
  pmer.set_unprojected_noise(dt);
  pmer.update_G();
  pmer.update_Hhat();
  pmer.compute_noise();
  pmer.compute_effective_kappa();
  pmer.compute_uc_forces();  
  pmer.add_external_force(pmer.dFdX,pmer.nuc_beads[0]);


  pmer.compute_tension({0.0,0.0,0.0});
  pmer.initial_integrate(dt);
  

  pmer.update_Hhat();
  pmer.compute_effective_kappa();
  pmer.compute_uc_forces();
  pmer.add_external_force(pmer.dFdX,pmer.nuc_beads[0]);

  
  pmer.compute_tension({0.0,0.0,0.0});
  int itermax = 20;
  
  int iterations = pmer.correct_tension(dt,pmer.atoms[0].R,itermax,1e-8);
  if (iterations > itermax) {
    std::cout << "too many iterations when correcting tension at time " << t
	      <<  ", retrying the step with new noise " << std::endl;
    for (int i = 0; i < pmer.get_Nbeads(); i++) 
      pmer.atoms[i].R = pmer.Rtmp[i];

    single_step(t,pmer,dt);
  }
  else {
    pmer.final_integrate(dt);
  
    pmer.compute_tangents_and_friction();
  }

  return;
}


std::string getLastLine(std::string filename);

void put_in_vectors(std::vector<std::vector<double>> & X_is,
		    std::string filename, MPI_Comm comm,int mpi_id);


void shift_pmer(std::vector<double> & X_i,
		SingleTetherWrap & pmer) {

  Eigen::Vector3d shift;

  int Nf = pmer.get_Nbeads() -1;

  shift(0) = pmer.atoms[Nf].R(0) - X_i[0];
  shift(1) = pmer.atoms[Nf].R(1) - X_i[1];
  shift(2) = pmer.atoms[Nf].R(2) - X_i[2];

  for (int i = 0; i <= Nf; i++) {
    pmer.atoms[i].R -= shift;
  }

  return;

}


void run(GlobalParams gp, psPDE::SolutionParams solparams,
	 const std::vector<std::string> &polymertypes,
	 const std::vector<std::vector<std::string>> & polymersplitvecs,
	 std::vector<std::vector<double>> &X_is) {


  // split the polymers between the processors.

  std::vector<std::unique_ptr<NoTetherWrap>> free_polys;
  std::vector<std::unique_ptr<SingleTetherWrap>> single_polys;
  std::vector<std::unique_ptr<DoubleTetherWrap>> double_polys;

  int stationary_nuc_count = X_is.size();
  
  
  // split up polymers of each type between all the processes, and assign them to a
  // nucleation site



  for (int i = 0; i < polymertypes.size(); i++) {

    if (i % gp.mpi_size == gp.id) {
      if (polymertypes.at(i) == "no_tether") {
	free_polys.push_back(
	 std::unique_ptr<NoTetherWrap>( new NoTetherWrap(polymersplitvecs.at(i),i,
							 polymertypes.at(i)
							 +std::to_string(i)))
			     );
      } else if (polymertypes.at(i) == "single_tether") {
	single_polys.push_back(
	 std::unique_ptr<SingleTetherWrap>( new SingleTetherWrap(polymersplitvecs.at(i),i,
								 polymertypes.at(i)
								 +std::to_string(i)))
			       );
      } else if (polymertypes.at(i) == "double_tether") {
	double_polys.push_back(
	 std::unique_ptr<DoubleTetherWrap>( new DoubleTetherWrap (polymersplitvecs.at(i),i,
								  polymertypes.at(i)
								  +std::to_string(i)))
			       );
      } else {
	throw std::runtime_error("Incompatible polymer type.");
      }
    }

  }


  // add dummy vectors to ensure the correct number of nucleation sites
  // in addition to the stationary sites
  for (int i = 0; i < polymertypes.size(); i++) {
    X_is.push_back({0,0,0});

  }


  psPDE::fftw_MPI_3Darray<double> phi(gp.comm,"concentration",gp.realspace);
  psPDE::fftw_MPI_3Darray<double> nonlinear(gp.comm,"chempotential",gp.realspace);


  psPDE::RandomPll rpll(gp.comm,gp.id,gp.seed,gp.mpi_size);
  
  psPDE::Integrator integrator(gp.comm,gp.fourier,rpll.get_processor_seed(),solparams,gp.dt);

  fftw_plan forward_phi, backward_phi;
  fftw_plan forward_nonlinear, backward_nonlinear;

  
  forward_phi = fftw_mpi_plan_dft_r2c_3d(gp.realspace.get_Nz(),gp.realspace.get_Ny(),
					 gp.realspace.get_Nx(),
					 phi.data(),
					 reinterpret_cast<fftw_complex*>
					 (integrator.ft_phi.data()),
					 gp.comm, FFTW_MPI_TRANSPOSED_OUT);
  
  backward_phi = fftw_mpi_plan_dft_c2r_3d(gp.realspace.get_Nz(),gp.realspace.get_Ny(),
					  gp.realspace.get_Nx(),
					  reinterpret_cast<fftw_complex*>
					  (integrator.ft_phi.data()),
					  phi.data(),gp.comm,FFTW_MPI_TRANSPOSED_IN);

  forward_nonlinear = fftw_mpi_plan_dft_r2c_3d(gp.realspace.get_Nz(),gp.realspace.get_Ny(),
					       gp.realspace.get_Nx(),
					       nonlinear.data(),
					       reinterpret_cast<fftw_complex*>
					       (integrator.ft_nonlinear.data()),
					       gp.comm, FFTW_MPI_TRANSPOSED_OUT);
  
  backward_nonlinear = fftw_mpi_plan_dft_c2r_3d(gp.realspace.get_Nz(),gp.realspace.get_Ny(),
						gp.realspace.get_Nx(),
						reinterpret_cast<fftw_complex*>
						(integrator.ft_nonlinear.data()),
						nonlinear.data(),gp.comm,
						FFTW_MPI_TRANSPOSED_IN);
  

  integrator.initialize(phi,gp.volFrac,gp.variance);

  
  std::vector<std::vector<double>> free_energy_derivs;

  std::string prefix = gp.solution_dump_file + std::string("_p") + std::to_string(gp.id) ;
  std::string fname_p = prefix + std::string("_") +  std::to_string(gp.startstep) +  std::string(".vti");

  std::string collection_name = prefix + std::string(".pvd");
  
  std::string complexprefix = gp.solution_dump_file + std::string("_complex")
    + std::string("_p") + std::to_string(gp.id) ;

  std::string complexfname_p = complexprefix + std::string("_") +  std::to_string(gp.startstep) +  std::string(".vti");

  std::string complexcollection_name = complexprefix + std::string(".pvd");
  
  psPDE::fftw_MPI_3Darray<double> modulus(gp.comm,integrator.ft_phi.get_name()+std::string("_mod"),
					  gp.fourier.get_positiveNx_grid());



  
  int running_average_count = 0;

  std::ofstream myfile;

  
  if (gp.restart_flag) {

    std::cout << "restarting!" << std::endl;
    // load in polymer data.

    for (auto &pmer : free_polys) {
      std::string poly_collection = gp.polymer_dump_file + pmer->name;
      std::string poly_fname = poly_collection + std::string("_") + std::to_string(gp.startstep)
	+ std::string(".vtp");

      BeadRodPmer::ioVTK::restartVTKcollection(poly_collection + std::string(".pvd"));
      BeadRodPmer::ioVTK::readVTKPolyData(*pmer,poly_fname);



      
      transfer_nucleation_site(X_is[pmer->nuc_index+stationary_nuc_count],
			       pmer->atoms[pmer->nuc_beads[0]].R,
			       gp.realspace.get_Lx(),gp.realspace.get_Ly(),
			       gp.realspace.get_Lz());
      
      pmer->compute_tangents_and_friction();
      
      pmer->set_Hhat();
      pmer->set_dCdlambda();
      pmer->set_G();

    }

    std::cout << "made it past free_polys! " << std::endl;
    
    
    for (auto &pmer : single_polys) {
      std::string poly_collection = gp.polymer_dump_file + pmer->name;
      std::string poly_fname = poly_collection + std::string("_") + std::to_string(gp.startstep)
	+ std::string(".vtp");

      BeadRodPmer::ioVTK::restartVTKcollection(poly_collection + std::string(".pvd"));
      BeadRodPmer::ioVTK::readVTKPolyData(*pmer,poly_fname);
      
      transfer_nucleation_site(X_is[pmer->nuc_index+stationary_nuc_count],
			       pmer->atoms[pmer->nuc_beads[0]].R,
			       gp.realspace.get_Lx(),gp.realspace.get_Ly(),
			       gp.realspace.get_Lz());
      

      pmer->compute_tangents_and_friction();
      
      pmer->set_Hhat();
      pmer->set_dCdlambda();
      pmer->set_G();
      
    }

    std::cout << "made it past single_polys! " << std::endl;
    
    for (auto &pmer : double_polys) {
      
      std::string poly_collection = gp.polymer_dump_file + pmer->name;
      std::string poly_fname = poly_collection + std::string("_") + std::to_string(gp.startstep)
	+ std::string(".vtp");


      BeadRodPmer::ioVTK::restartVTKcollection(poly_collection + std::string(".pvd"));
      BeadRodPmer::ioVTK::readVTKPolyData(*pmer,poly_fname);

      transfer_nucleation_site(X_is[pmer->nuc_index+stationary_nuc_count],
			       pmer->atoms[pmer->nuc_beads[0]].R,
			       gp.realspace.get_Lx(),gp.realspace.get_Ly(),
			       gp.realspace.get_Lz());
      

      pmer->compute_tangents_and_friction();
      
      pmer->set_Hhat();
      pmer->set_dCdlambda();
      pmer->set_G();
      
    }

    std::cout << "made it past double_polys! " << std::endl;
    
    for (int i = 0; i < polymertypes.size(); i++)
      MPI_Bcast(X_is[i+stationary_nuc_count].data(),3,MPI_DOUBLE,i % gp.mpi_size,gp.comm);


    std::cout << "made it past the broadcasting of all nucleation points! " << std::endl;

    psPDE::ioVTK::restartVTKcollection(collection_name);
    std::cout << "made it past the restart vtk real collection! " << std::endl;
    
    psPDE::ioVTK::restartVTKcollection(complexcollection_name);
    std::cout << "made it past the restart vtk complex collection! " << std::endl;
    
    
    psPDE::ioVTK::readVTKImageData({&phi},fname_p);
    std::cout << "made it past the restart vtk read image data! " << std::endl;
    
    
    for (int i = 0; i < integrator.ft_phi.axis_size(0); i++) {
      for (int j = 0; j < integrator.ft_phi.axis_size(1); j++) {
	for (int k = 0; k < integrator.ft_phi.axis_size(2); k++) {
	  modulus(i,j,k) = 0.0;
	}
      }
    }


    std::cout << "made it past solution initialise! " << std::endl;
    // open thermo data so we can append to it

    if (gp.id == 0) {
      myfile.open(gp.thermo_file,std::ios::app);
      
    }

    std::cout << "made it past appending data! " << std::endl;
    
  } else {

    for (auto &pmer : free_polys) {
      
      std::string poly_collection = gp.polymer_dump_file + pmer->name;
      
      BeadRodPmer::ioVTK::writeVTKcollectionHeader(poly_collection + std::string(".pvd"));

      transfer_nucleation_site(X_is[pmer->nuc_index+stationary_nuc_count],
			       pmer->atoms[pmer->nuc_beads[0]].R,
			       gp.realspace.get_Lx(),gp.realspace.get_Ly(),
			       gp.realspace.get_Lz());
      
      pmer->compute_tangents_and_friction();
      
      pmer->set_Hhat();
      pmer->set_dCdlambda();
      pmer->set_G();
      
    }
    
    for (auto &pmer : single_polys) {

      std::string poly_collection = gp.polymer_dump_file + pmer->name;
      
      BeadRodPmer::ioVTK::writeVTKcollectionHeader(poly_collection + std::string(".pvd"));
      
      transfer_nucleation_site(X_is[pmer->nuc_index+stationary_nuc_count],
			       pmer->atoms[pmer->nuc_beads[0]].R,
			       gp.realspace.get_Lx(),gp.realspace.get_Ly(),
			       gp.realspace.get_Lz());
      
      
      pmer->compute_tangents_and_friction();
      
      pmer->set_Hhat();
      pmer->set_dCdlambda();
      pmer->set_G();
    }
    
    for (auto &pmer : double_polys) {

      std::string poly_collection = gp.polymer_dump_file + pmer->name;
      
      BeadRodPmer::ioVTK::writeVTKcollectionHeader(poly_collection + std::string(".pvd"));
      
      transfer_nucleation_site(X_is[pmer->nuc_index+stationary_nuc_count],
			       pmer->atoms[pmer->nuc_beads[0]].R,
			       gp.realspace.get_Lx(),gp.realspace.get_Ly(),
			       gp.realspace.get_Lz());
      
      pmer->compute_tangents_and_friction();
      
      pmer->set_Hhat();
      pmer->set_dCdlambda();
      pmer->set_G();
      
    }
    for (int i = 0; i < polymertypes.size(); i++)
      MPI_Bcast(X_is[i+stationary_nuc_count].data(),3,MPI_DOUBLE,i % gp.mpi_size,gp.comm);
        


    // equilibrate polymers before nucleating
    
    double teq = 0;
    for (int it = 1; it <= gp.polymer_equilibration; it++) {

      teq += gp.dt;
      // update the polymers
      for (auto &pmer : free_polys)  {

	single_step(teq,*pmer,gp.dt);
	transfer_nucleation_site(X_is[pmer->nuc_index+stationary_nuc_count],
				 pmer->atoms[pmer->nuc_beads[0]].R,
				 gp.realspace.get_Lx(),gp.realspace.get_Ly(),
				 gp.realspace.get_Lz());
	
      }
      for (auto &pmer : single_polys) {

	single_step(teq,*pmer,gp.dt);
	transfer_nucleation_site(X_is[pmer->nuc_index+stationary_nuc_count],
				 pmer->atoms[pmer->nuc_beads[0]].R,
				 gp.realspace.get_Lx(),gp.realspace.get_Ly(),
				 gp.realspace.get_Lz());
	
      }
      for (auto &pmer : double_polys) {

	single_step(teq,*pmer,gp.dt);
	transfer_nucleation_site(X_is[pmer->nuc_index+stationary_nuc_count],
				 pmer->atoms[pmer->nuc_beads[0]].R,
				 gp.realspace.get_Lx(),gp.realspace.get_Ly(),
				 gp.realspace.get_Lz());
	
      }

      // and also their nucleation site positions.
      for (int i = 0; i < polymertypes.size(); i++)
	MPI_Bcast(X_is[i+stationary_nuc_count].data(),3,MPI_DOUBLE,i % gp.mpi_size,gp.comm);
      
    }
    


    
    // write equilibrated state.
    for (auto &pmer : free_polys) {
      
      std::string poly_collection = gp.polymer_dump_file + pmer->name;
      std::string poly_fname = poly_collection + std::string("_") + std::to_string(0)
	+ std::string(".vtp");
      BeadRodPmer::ioVTK::writeVTKPolyData(poly_fname,*pmer);
      BeadRodPmer::ioVTK::writeVTKcollectionMiddle(poly_collection+ std::string(".pvd")
						   ,poly_fname,gp.starttime);
    }


    // write equilibrated state.
    for (auto &pmer : single_polys) {
      
      std::string poly_collection = gp.polymer_dump_file + pmer->name;
      std::string poly_fname = poly_collection + std::string("_") + std::to_string(0)
	+ std::string(".vtp");
      BeadRodPmer::ioVTK::writeVTKPolyData(poly_fname,*pmer);
      BeadRodPmer::ioVTK::writeVTKcollectionMiddle(poly_collection+ std::string(".pvd")
						   ,poly_fname,gp.starttime);
    }

    // write equilibrated state.
    for (auto &pmer : double_polys) {
      
      std::string poly_collection = gp.polymer_dump_file + pmer->name;
      std::string poly_fname = poly_collection + std::string("_") + std::to_string(0)
	+ std::string(".vtp");
      BeadRodPmer::ioVTK::writeVTKPolyData(poly_fname,*pmer);
      BeadRodPmer::ioVTK::writeVTKcollectionMiddle(poly_collection+ std::string(".pvd")
						   ,poly_fname,gp.starttime);
    }



    psPDE::ioVTK::writeVTKcollectionHeader(collection_name);
    psPDE::ioVTK::writeVTKcollectionHeader(complexcollection_name);
    
    fftw_execute(forward_phi);
    integrator.ft_phi.mod(modulus);
    
    double norm = 1.0/(integrator.ft_phi.grid.get_Nx()*integrator.ft_phi.grid.get_Ny()
		       *integrator.ft_phi.grid.get_Nz());
    
    for (int i = 0; i < integrator.ft_phi.axis_size(0); i++) {
      for (int j = 0; j < integrator.ft_phi.axis_size(1); j++) {
	for (int k = 0; k < integrator.ft_phi.axis_size(2); k++) {
	  integrator.ft_phi(i,j,k) = integrator.ft_phi(i,j,k)*norm;
	  modulus(i,j,k) = modulus(i,j,k)*norm;
	}
      }
    }
    if (gp.id == 0) {
      modulus(0,0,0) = 0.0;
    }

    fftw_execute(backward_phi);
    
    psPDE::ioVTK::writeVTKImageData(fname_p,{&phi},gp.realspace);
    psPDE::ioVTK::writeVTKImageData(complexfname_p,{&modulus},modulus.grid);
    
    psPDE::ioVTK::writeVTKcollectionMiddle(collection_name,fname_p,gp.starttime);
    psPDE::ioVTK::writeVTKcollectionMiddle(complexcollection_name,complexfname_p,gp.starttime);

    
    
    if (gp.id == 0) {
      myfile.open(gp.thermo_file);
      myfile << "# t ";
      for (unsigned index = 0; index < X_is.size() ; index ++) {
	myfile << "\t (X_" << index << ")_x " << "(X_" << index << ")_y"
	       << "(X_" << index << ")_z";
      }
      myfile << "\t F(X) ";
      for (unsigned index = 0; index < X_is.size() ; index ++) {
	myfile << "\t (dF/dX_" << index << ")_x " << "(dF/dX_" << index << ")_y"
	       << "(dF/dX_" << index << ")_z";
      }
      myfile << std::endl;
      
    }
    
    
  }

  std::cout << " polymers and solution are initialized and/or equilibrated." << std::endl;
  
  // main loop!
  
  psPDE::TimeStep timestep(gp.comm,gp.mpi_size,gp.id,integrator.ft_phi.axis_size(0),
		    integrator.ft_phi.axis_size(1));

  double free_energy;

  double t = gp.starttime;


  for (int it = 1+gp.startstep; it <= gp.steps+gp.startstep; it ++) {

    if (abs(t-1.28463) < 1e-12 || abs(t-1.284625) < 1e-12 ) {

      std::cout << "going to get stuck on this step. " << std::endl;
      // compute nl(t) given phi(t), and also free energy derivatives for the different
      //   nucleation sites.
      free_energy_derivs = integrator.nonlinear(nonlinear,phi,X_is,free_energy);
      std::cout << " did the integration step. " << std::endl;

      if (gp.id == 0 ) {
	std::cout.precision(10);
	std::cout << std::fixed;
	std::cout << "printing here " << t << std::endl;
      }
      
      if (gp.id == 0 && it % gp.thermo_every == 0) {
	myfile << t;
	
	for (unsigned index = 0; index < X_is.size() ; index ++) {
	  
	  myfile << "\t " << X_is.at(index).at(0) << "\t "
		 << X_is.at(index).at(1) << "\t "
		 << X_is.at(index).at(2);
	}
	
	myfile << "\t " << free_energy;
	for (unsigned index = 0; index < free_energy_derivs.size() ; index ++) {
	  
	  myfile << "\t " << free_energy_derivs.at(index).at(0) << "\t "
		 << free_energy_derivs.at(index).at(1) << "\t "
		 << free_energy_derivs.at(index).at(2);
	}
	
	
	myfile << std::endl;
	
      }
      
      t += gp.dt;
      std::cout << " saved the free energy. " << std::endl;
      // update the nucleation sites via the polymers.
      
      for (auto &pmer : free_polys)  {
	transfer_free_energy(pmer->dFdX,free_energy_derivs,
			     pmer->nuc_index+stationary_nuc_count);
	single_step(t,*pmer,gp.dt);
	transfer_nucleation_site(X_is[pmer->nuc_index+stationary_nuc_count],
				 pmer->atoms[pmer->nuc_beads[0]].R,
				 gp.realspace.get_Lx(),gp.realspace.get_Ly(),
				 gp.realspace.get_Lz());
	
	
	if ( it % gp.polymer_dump_every == 0) {
	  std::string poly_collection = gp.polymer_dump_file + pmer->name;
	  std::string poly_fname = poly_collection + std::string("_") + std::to_string(it)
	    + std::string(".vtp");
	  BeadRodPmer::ioVTK::writeVTKPolyData(poly_fname,*pmer);
	  BeadRodPmer::ioVTK::writeVTKcollectionMiddle(poly_collection+ std::string(".pvd"),
						       poly_fname,t);
	}
	
      }
      
      
      for (auto &pmer : single_polys)  {

	std::string poly_collection = gp.polymer_dump_file + pmer->name;
	std::string poly_fname = poly_collection + std::string("_") + std::to_string(it)
	  + std::string(".vtp");
	BeadRodPmer::ioVTK::writeVTKPolyData(poly_fname,*pmer);
	BeadRodPmer::ioVTK::writeVTKcollectionMiddle(poly_collection+ std::string(".pvd"),
						     poly_fname,t);
	
	
	transfer_free_energy(pmer->dFdX,free_energy_derivs,
			     pmer->nuc_index+stationary_nuc_count);

	single_step(t,*pmer,gp.dt);
	transfer_nucleation_site(X_is[pmer->nuc_index+stationary_nuc_count],
				 pmer->atoms[pmer->nuc_beads[0]].R,
				 gp.realspace.get_Lx(),gp.realspace.get_Ly(),
				 gp.realspace.get_Lz());
	
	
	
      }

      std::cout << " got past the single_polys on id " << gp.id << std::endl;                  
      for (auto &pmer : double_polys)  {
	transfer_free_energy(pmer->dFdX,free_energy_derivs,
			     pmer->nuc_index+stationary_nuc_count);
	single_step(t,*pmer,gp.dt);
	transfer_nucleation_site(X_is[pmer->nuc_index+stationary_nuc_count],
				 pmer->atoms[pmer->nuc_beads[0]].R,
				 gp.realspace.get_Lx(),gp.realspace.get_Ly(),
				 gp.realspace.get_Lz());
	
	
	if ( it % gp.polymer_dump_every == 0) {
	  std::string poly_collection = gp.polymer_dump_file + pmer->name;
	  std::string poly_fname = poly_collection + std::string("_") + std::to_string(it)
	    + std::string(".vtp");
	  BeadRodPmer::ioVTK::writeVTKPolyData(poly_fname,*pmer);
	  BeadRodPmer::ioVTK::writeVTKcollectionMiddle(poly_collection+ std::string(".pvd"),
						       poly_fname,t);
	}
	
      }
      
      std::cout << " got past the double_polys. " << std::endl;      
      
      // update the volume fraction phi
      
      fftw_execute(forward_phi);
      fftw_execute(forward_nonlinear);
      std::cout << "survived the forward fft. " << std::endl;      
      
      timestep.update(t,integrator);

      std::cout << "survived timestep update. " << std::endl;
      
      integrator.ft_phi.running_mod(modulus);
      running_average_count += 1;
      
      if (it % gp.solution_dump_every == 0) {
	
	if (gp.id == 0) {
	  std::cout << " saving at t = " << t << std::endl;
	}
	
	
	fname_p = prefix + std::string("_") +  std::to_string(it) +  std::string(".vti");
	complexfname_p = complexprefix + std::string("_") +  std::to_string(it) +  std::string(".vti");
	modulus /= running_average_count;
	if (gp.id == 0) {
	  modulus(0,0,0) = 0.0;
	}
	
	running_average_count = 0;
	
	psPDE::ioVTK::writeVTKImageData(complexfname_p,{&modulus},modulus.grid);
	psPDE::ioVTK::writeVTKcollectionMiddle(complexcollection_name,complexfname_p,t);
	fftw_execute(backward_phi); // get phi(t+dt)
	
	integrator.initialize(modulus,0,0);
	
	psPDE::ioVTK::writeVTKImageData(fname_p,{&phi},phi.grid);
	psPDE::ioVTK::writeVTKcollectionMiddle(collection_name,fname_p,t);
      } else {
	fftw_execute(backward_phi); // get phi(t+dt)
      }
      
      // broadcast the new nucleation sites across the processes
      for (int i = 0; i < polymertypes.size(); i++)
	MPI_Bcast(X_is[i+stationary_nuc_count].data(),3,MPI_DOUBLE,i % gp.mpi_size,gp.comm);
      
      
      
    } else {
      
      // compute nl(t) given phi(t), and also free energy derivatives for the different
      //   nucleation sites.
      free_energy_derivs = integrator.nonlinear(nonlinear,phi,X_is,free_energy); 
      
      
      if (gp.id == 0 && it % gp.thermo_every == 0) {
	myfile << t;
	
	for (unsigned index = 0; index < X_is.size() ; index ++) {
	  
	  myfile << "\t " << X_is.at(index).at(0) << "\t "
		 << X_is.at(index).at(1) << "\t "
		 << X_is.at(index).at(2);
	}
	
	myfile << "\t " << free_energy;
	for (unsigned index = 0; index < free_energy_derivs.size() ; index ++) {
	  
	  myfile << "\t " << free_energy_derivs.at(index).at(0) << "\t "
		 << free_energy_derivs.at(index).at(1) << "\t "
		 << free_energy_derivs.at(index).at(2);
	}
	
	
	myfile << std::endl;
	
      }
      
      t += gp.dt;
      
      // update the nucleation sites via the polymers.
      
      for (auto &pmer : free_polys)  {
	transfer_free_energy(pmer->dFdX,free_energy_derivs,
			     pmer->nuc_index+stationary_nuc_count);
	single_step(t,*pmer,gp.dt);
	transfer_nucleation_site(X_is[pmer->nuc_index+stationary_nuc_count],
				 pmer->atoms[pmer->nuc_beads[0]].R,
				 gp.realspace.get_Lx(),gp.realspace.get_Ly(),
				 gp.realspace.get_Lz());
	
	
	if ( it % gp.polymer_dump_every == 0) {
	  std::string poly_collection = gp.polymer_dump_file + pmer->name;
	  std::string poly_fname = poly_collection + std::string("_") + std::to_string(it)
	    + std::string(".vtp");
	  BeadRodPmer::ioVTK::writeVTKPolyData(poly_fname,*pmer);
	  BeadRodPmer::ioVTK::writeVTKcollectionMiddle(poly_collection+ std::string(".pvd"),
						       poly_fname,t);
	}
	
      }
      
      
      for (auto &pmer : single_polys)  {
	transfer_free_energy(pmer->dFdX,free_energy_derivs,
			     pmer->nuc_index+stationary_nuc_count);
	single_step(t,*pmer,gp.dt);
	transfer_nucleation_site(X_is[pmer->nuc_index+stationary_nuc_count],
				 pmer->atoms[pmer->nuc_beads[0]].R,
				 gp.realspace.get_Lx(),gp.realspace.get_Ly(),
				 gp.realspace.get_Lz());
	
	
	if ( it % gp.polymer_dump_every == 0) {
	  
	  
	  std::string poly_collection = gp.polymer_dump_file + pmer->name;
	  std::string poly_fname = poly_collection + std::string("_") + std::to_string(it)
	    + std::string(".vtp");
	  BeadRodPmer::ioVTK::writeVTKPolyData(poly_fname,*pmer);
	  BeadRodPmer::ioVTK::writeVTKcollectionMiddle(poly_collection+ std::string(".pvd"),
						       poly_fname,t);
	}
	
      }
      
      
      for (auto &pmer : double_polys)  {
	transfer_free_energy(pmer->dFdX,free_energy_derivs,
			     pmer->nuc_index+stationary_nuc_count);
	single_step(t,*pmer,gp.dt);
	transfer_nucleation_site(X_is[pmer->nuc_index+stationary_nuc_count],
				 pmer->atoms[pmer->nuc_beads[0]].R,
				 gp.realspace.get_Lx(),gp.realspace.get_Ly(),
				 gp.realspace.get_Lz());
	
	
	if ( it % gp.polymer_dump_every == 0) {
	  std::string poly_collection = gp.polymer_dump_file + pmer->name;
	  std::string poly_fname = poly_collection + std::string("_") + std::to_string(it)
	    + std::string(".vtp");
	  BeadRodPmer::ioVTK::writeVTKPolyData(poly_fname,*pmer);
	  BeadRodPmer::ioVTK::writeVTKcollectionMiddle(poly_collection+ std::string(".pvd"),
						       poly_fname,t);
	}
	
      }
      
      
      
      // update the volume fraction phi
      
      fftw_execute(forward_phi);
      fftw_execute(forward_nonlinear);
      
      
      timestep.update(t,integrator);
      
      integrator.ft_phi.running_mod(modulus);
      running_average_count += 1;
      
      if (it % gp.solution_dump_every == 0) {
	
	if (gp.id == 0) {
	  std::cout << " saving at t = " << t << std::endl;
	}
	
	
	fname_p = prefix + std::string("_") +  std::to_string(it) +  std::string(".vti");
	complexfname_p = complexprefix + std::string("_") +  std::to_string(it) +  std::string(".vti");
	modulus /= running_average_count;
	if (gp.id == 0) {
	  modulus(0,0,0) = 0.0;
	}
	
	running_average_count = 0;
	
	psPDE::ioVTK::writeVTKImageData(complexfname_p,{&modulus},modulus.grid);
	psPDE::ioVTK::writeVTKcollectionMiddle(complexcollection_name,complexfname_p,t);
	fftw_execute(backward_phi); // get phi(t+dt)
	
	integrator.initialize(modulus,0,0);
	
	psPDE::ioVTK::writeVTKImageData(fname_p,{&phi},phi.grid);
	psPDE::ioVTK::writeVTKcollectionMiddle(collection_name,fname_p,t);
      } else {
	fftw_execute(backward_phi); // get phi(t+dt)
      }
      
      // broadcast the new nucleation sites across the processes
      for (int i = 0; i < polymertypes.size(); i++)
	MPI_Bcast(X_is[i+stationary_nuc_count].data(),3,MPI_DOUBLE,i % gp.mpi_size,gp.comm);
      
      
      
    }

  }
    
  if (gp.id == 0) {
    myfile.close();
  }
  
  psPDE::ioVTK::writeVTKcollectionFooter(collection_name);
  psPDE::ioVTK::writeVTKcollectionFooter(complexcollection_name);
  
  
  for (auto &pmer : free_polys)  {
    
    std::string poly_collection = gp.polymer_dump_file + pmer->name + std::string(".pvd");
    
    BeadRodPmer::ioVTK::writeVTKcollectionFooter(poly_collection);
    
  }
  
  for (auto &pmer : single_polys)  {
    
    std::string poly_collection = gp.polymer_dump_file + pmer->name + std::string(".pvd");
    
    BeadRodPmer::ioVTK::writeVTKcollectionFooter(poly_collection);
    
  }
  
  
  for (auto &pmer : double_polys)  {
    
    std::string poly_collection = gp.polymer_dump_file + pmer->name + std::string(".pvd");
    
    BeadRodPmer::ioVTK::writeVTKcollectionFooter(poly_collection);
    
  }
  
  
  
  fftw_destroy_plan(forward_phi);
  fftw_destroy_plan(backward_phi);
  fftw_destroy_plan(forward_nonlinear);
  fftw_destroy_plan(backward_nonlinear);
  
  
  return;
}



std::string getLastLine(std::string filename)
{

  std::ifstream myfile (filename);
  std::string lastline;

  if (myfile) {

    // assumes that last line is just a newline char, hence the -2 below (vs -1)
    myfile.seekg(-2,myfile.end);

    char ch;

    myfile.get(ch);

    while (ch != '\n' && myfile.tellg() >1) {
    
      myfile.seekg(-2,myfile.cur);
      myfile.get(ch);
    }




    std::getline(myfile,lastline);

  } else {
    throw std::runtime_error("File " + std::string(filename) + " does not exist");
  }

  myfile.close();
  return lastline;
}


void put_in_vectors(std::vector<std::vector<double>> & X_is,
		    std::string filename, MPI_Comm comm,int mpi_id)
{

  std::string finalline;
  int signal = 0;
  if (mpi_id == 0) {

    try {
      finalline = getLastLine(filename);
    }
    catch (const std::runtime_error & error) {
      signal = -1;
    }

  }

  MPI_Bcast(&signal,1,MPI_INT,0,comm);

  if (signal == -1) {
    throw std::runtime_error("File " + std::string(filename) + " does not exist");    
  }

  std::istringstream iss(finalline);
  
  std::string subs;
  if (mpi_id == 0) {
    iss >> subs;
  }

  for (auto &cmp : X_is) {
    for (int i = 0; i < 3; i++) {
      if (mpi_id == 0) {
	iss >> cmp.at(i);
      }

      MPI_Bcast(&cmp.at(i),1,MPI_DOUBLE,0,comm);
    }
  }

  return;
  
}
