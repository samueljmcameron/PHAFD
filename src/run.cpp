#include <fstream>
#include <memory>
#include <chrono>

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

std::string getLastLine(std::string filename);

void put_in_vectors(std::vector<std::vector<double>> & X_is,
		    std::string filename, MPI_Comm comm,int mpi_id);






void split_polymers_across_processors(const std::vector<std::string>
				      &polymertypes,
				      const std::vector<std::vector<std::string>>
				      & polymersplitvecs,
				      std::vector<std::vector<double>> &X_is,
				      std::vector<int> & beads_per_pmer,
				      std::vector<std::unique_ptr<NoTetherWrap>>
				      &free_polys,
				      std::vector<std::unique_ptr<SingleTetherWrap>> 
				      &single_polys,
				      std::vector<std::unique_ptr<DoubleTetherWrap>>
				      &double_polys,int mpi_size,int id)
{
  
  int running_total_beads = 0; 

  for (std::vector<std::string>::size_type i = 0; i < polymertypes.size(); i++) {

    if (polymertypes.at(i) == "no_tether") {
      
      std::unique_ptr<NoTetherWrap>
	tmpobj(new NoTetherWrap(polymersplitvecs.at(i),i,
				running_total_beads,
				polymertypes.at(i)+std::to_string(i)));
      
      beads_per_pmer.push_back(tmpobj->nuc_beads.size());

      if (i % mpi_size == id) {
	free_polys.push_back( std::move(tmpobj));
	free_polys.back()->dFdX_is.resize(beads_per_pmer.back());
	for (auto &dFdX : free_polys.back()->dFdX_is)
	  dFdX = {0,0,0};
	  
      }

    } else if (polymertypes.at(i) == "single_tether") {

      std::unique_ptr<SingleTetherWrap>
	tmpobj(new SingleTetherWrap(polymersplitvecs.at(i),i,
				    running_total_beads,
				    polymertypes.at(i)+std::to_string(i)));
      
      beads_per_pmer.push_back(tmpobj->nuc_beads.size());

      
      if (i % mpi_size == id) {
	single_polys.push_back(std::move(tmpobj) );
	single_polys.back()->dFdX_is.resize(beads_per_pmer.back());
	for (auto &dFdX : single_polys.back()->dFdX_is)
	  dFdX = {0,0,0};
	
      }

    } else if (polymertypes.at(i) == "double_tether") {

      std::unique_ptr<DoubleTetherWrap>
	tmpobj(new DoubleTetherWrap(polymersplitvecs.at(i),i,
				    running_total_beads,
				    polymertypes.at(i)+std::to_string(i)));
      
      beads_per_pmer.push_back(tmpobj->nuc_beads.size());
      
      if (i % mpi_size == id) {
	double_polys.push_back(std::move(tmpobj));
	double_polys.back()->dFdX_is.resize(beads_per_pmer.back());
	for (auto &dFdX : double_polys.back()->dFdX_is)
	  dFdX = {0,0,0};
      }

    } else {
      throw std::runtime_error("Incompatible polymer type.");
    }
    running_total_beads += beads_per_pmer.back();
  }

  
  for (int index = 0; index < running_total_beads; index ++) 
    X_is.push_back({0,0,0});
      


  return;
}

void broadcast_X_is(int stationary_nuc_count, const std::vector<int> &beads_per_pmer,
		    std::vector<std::vector<double>> &X_is,int mpi_size,
		    MPI_Comm comm)
{
  int startpoint = stationary_nuc_count;    
  for (int index = 0; index <  beads_per_pmer.size(); index ++) {
    
    if (index != 0) {
      startpoint += beads_per_pmer.at(index-1);
    }
    int ni = beads_per_pmer.at(index);
    
    for (int ipos = 0; ipos < ni; ipos ++)
      MPI_Bcast(X_is[startpoint + ipos].data(),3,MPI_DOUBLE,
		index % mpi_size,comm);
  }
  return;
}

void run(GlobalParams gp, psPDE::SolutionParams solparams,
	 const std::vector<std::string> &polymertypes,
	 const std::vector<std::vector<std::string>> & polymersplitvecs,
	 std::vector<std::vector<double>> &X_is) {

  
  std::vector<std::unique_ptr<NoTetherWrap>> free_polys;
  std::vector<std::unique_ptr<SingleTetherWrap>> single_polys;
  std::vector<std::unique_ptr<DoubleTetherWrap>> double_polys;

  int stationary_nuc_count = X_is.size();
  
  std::vector<int> beads_per_pmer;
  
  // split the polymers between the processors.
  split_polymers_across_processors(polymertypes,polymersplitvecs,X_is,
				   beads_per_pmer,free_polys,single_polys,
				   double_polys,gp.mpi_size,gp.id);

  
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
    if (gp.id == 0) {
      std::cout << "restarting!" << std::endl;
    }
    // load in polymer data.

    for (auto &pmer : free_polys) {
      std::string poly_collection = gp.polymer_dump_file + pmer->name;
      std::string poly_fname = poly_collection + std::string("_") + std::to_string(gp.startstep)
	+ std::string(".vtp");

      BeadRodPmer::ioVTK::restartVTKcollection(poly_collection + std::string(".pvd"));
      BeadRodPmer::ioVTK::readVTKPolyData(*pmer,poly_fname);


      for (int index = 0; index < beads_per_pmer.at(pmer->number); index ++) {
	transfer_nucleation_site(X_is[pmer->nuc_index+index+stationary_nuc_count],
				 pmer->atoms[pmer->nuc_beads[index]].R,
				 gp.realspace.get_Lx(),gp.realspace.get_Ly(),
				 gp.realspace.get_Lz());
      }
      
      pmer->compute_tangents_and_friction();
      
      pmer->set_Hhat();
      pmer->set_dCdlambda();
      pmer->set_G();

    }


    
    
    for (auto &pmer : single_polys) {
      std::string poly_collection = gp.polymer_dump_file + pmer->name;
      std::string poly_fname = poly_collection + std::string("_") + std::to_string(gp.startstep)
	+ std::string(".vtp");

      BeadRodPmer::ioVTK::restartVTKcollection(poly_collection + std::string(".pvd"));
      BeadRodPmer::ioVTK::readVTKPolyData(*pmer,poly_fname);

      for (int index = 0; index < beads_per_pmer.at(pmer->number); index ++) {
	transfer_nucleation_site(X_is[pmer->nuc_index+index+stationary_nuc_count],
				 pmer->atoms[pmer->nuc_beads[index]].R,
				 gp.realspace.get_Lx(),gp.realspace.get_Ly(),
				 gp.realspace.get_Lz());
      }



      pmer->compute_tangents_and_friction();
      
      pmer->set_Hhat();
      pmer->set_dCdlambda();
      pmer->set_G();
      
    }


    
    for (auto &pmer : double_polys) {
      
      std::string poly_collection = gp.polymer_dump_file + pmer->name;
      std::string poly_fname = poly_collection + std::string("_") + std::to_string(gp.startstep)
	+ std::string(".vtp");


      BeadRodPmer::ioVTK::restartVTKcollection(poly_collection + std::string(".pvd"));
      BeadRodPmer::ioVTK::readVTKPolyData(*pmer,poly_fname);


      for (int index = 0; index < beads_per_pmer.at(pmer->number); index ++) {
	transfer_nucleation_site(X_is[pmer->nuc_index+index+stationary_nuc_count],
				 pmer->atoms[pmer->nuc_beads[index]].R,
				 gp.realspace.get_Lx(),gp.realspace.get_Ly(),
				 gp.realspace.get_Lz());
      }

      pmer->compute_tangents_and_friction();
      
      pmer->set_Hhat();
      pmer->set_dCdlambda();
      pmer->set_G();
      
    }
    

    broadcast_X_is(stationary_nuc_count,beads_per_pmer,X_is,gp.mpi_size,gp.comm);

    psPDE::ioVTK::restartVTKcollection(collection_name,gp.comm);

    
    psPDE::ioVTK::restartVTKcollection(complexcollection_name,gp.comm);

    
    
    psPDE::ioVTK::readVTKImageData({&phi},fname_p);

    
    
    for (int i = 0; i < integrator.ft_phi.axis_size(0); i++) {
      for (int j = 0; j < integrator.ft_phi.axis_size(1); j++) {
	for (int k = 0; k < integrator.ft_phi.axis_size(2); k++) {
	  modulus(i,j,k) = 0.0;
	}
      }
    }



    // open thermo data so we can append to it

    if (gp.id == 0) {
      myfile.open(gp.thermo_file,std::ios::app);
      
    }


    
  } else {

    for (auto &pmer : free_polys) {
      
      std::string poly_collection = gp.polymer_dump_file + pmer->name;
      
      BeadRodPmer::ioVTK::writeVTKcollectionHeader(poly_collection + std::string(".pvd"));

      for (int index = 0; index < beads_per_pmer.at(pmer->number); index ++) {
	transfer_nucleation_site(X_is[pmer->nuc_index+index+stationary_nuc_count],
				 pmer->atoms[pmer->nuc_beads[index]].R,
				 gp.realspace.get_Lx(),gp.realspace.get_Ly(),
				 gp.realspace.get_Lz());
      }
      
      pmer->compute_tangents_and_friction();
      
      pmer->set_Hhat();
      pmer->set_dCdlambda();
      pmer->set_G();
      
    }
    
    for (auto &pmer : single_polys) {

      std::string poly_collection = gp.polymer_dump_file + pmer->name;
      
      BeadRodPmer::ioVTK::writeVTKcollectionHeader(poly_collection + std::string(".pvd"));
      
      for (int index = 0; index < beads_per_pmer.at(pmer->number); index ++) {
	transfer_nucleation_site(X_is[pmer->nuc_index+index+stationary_nuc_count],
				 pmer->atoms[pmer->nuc_beads[index]].R,
				 gp.realspace.get_Lx(),gp.realspace.get_Ly(),
				 gp.realspace.get_Lz());
      }
      
      pmer->compute_tangents_and_friction();
      
      pmer->set_Hhat();
      pmer->set_dCdlambda();
      pmer->set_G();
    }
    
    for (auto &pmer : double_polys) {

      std::string poly_collection = gp.polymer_dump_file + pmer->name;
      
      BeadRodPmer::ioVTK::writeVTKcollectionHeader(poly_collection + std::string(".pvd"));
      
      for (int index = 0; index < beads_per_pmer.at(pmer->number); index ++) {
	transfer_nucleation_site(X_is[pmer->nuc_index+index+stationary_nuc_count],
				 pmer->atoms[pmer->nuc_beads[index]].R,
				 gp.realspace.get_Lx(),gp.realspace.get_Ly(),
				 gp.realspace.get_Lz());
      }
      
      pmer->compute_tangents_and_friction();
      
      pmer->set_Hhat();
      pmer->set_dCdlambda();
      pmer->set_G();
      
    }

    broadcast_X_is(stationary_nuc_count,beads_per_pmer,X_is,gp.mpi_size,gp.comm);

    // equilibrate polymers before nucleating
    
    double teq = 0;
    for (int it = 1; it <= gp.polymer_equilibration; it++) {

      teq += gp.dt;
      // update the polymers
      for (auto &pmer : free_polys)  {

	pmer->single_step(teq,gp.dt,pmer->dFdX_is);

	for (int index = 0; index < beads_per_pmer.at(pmer->number); index ++) {
	  transfer_nucleation_site(X_is[pmer->nuc_index+index+stationary_nuc_count],
				   pmer->atoms[pmer->nuc_beads[index]].R,
				   gp.realspace.get_Lx(),gp.realspace.get_Ly(),
				   gp.realspace.get_Lz());
	}
	

      }
      for (auto &pmer : single_polys) {

	pmer->single_step(teq,gp.dt,pmer->dFdX_is);	
	
	for (int index = 0; index < beads_per_pmer.at(pmer->number); index ++) {
	  transfer_nucleation_site(X_is[pmer->nuc_index+index+stationary_nuc_count],
				   pmer->atoms[pmer->nuc_beads[index]].R,
				   gp.realspace.get_Lx(),gp.realspace.get_Ly(),
				   gp.realspace.get_Lz());
	}

      }
      for (auto &pmer : double_polys) {
	
	pmer->single_step(teq,gp.dt,pmer->dFdX_is);	

	for (int index = 0; index < beads_per_pmer.at(pmer->number); index ++) {
	  transfer_nucleation_site(X_is[pmer->nuc_index+index+stationary_nuc_count],
				   pmer->atoms[pmer->nuc_beads[index]].R,
				   gp.realspace.get_Lx(),gp.realspace.get_Ly(),
				   gp.realspace.get_Lz());
	}

	
      }

      // and also their nucleation site positions.
      broadcast_X_is(stationary_nuc_count,beads_per_pmer,X_is,gp.mpi_size,gp.comm);
      
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

  std::cout << " all polymers and solution are initialized and/or equilibrated on process"
	    << gp.id << "." << std::endl;
  
  // main loop!
  
  psPDE::TimeStep timestep(gp.comm,gp.mpi_size,gp.id,integrator.ft_phi.axis_size(0),
			   integrator.ft_phi.axis_size(1));

  double free_energy;

  double t = gp.starttime;


  using Clock = std::chrono::steady_clock;
  using namespace std::literals;
  auto constexpr chronoitvl = 1.0s/60.0;
  using duration = std::chrono::duration<double>;
  using time_point = std::chrono::time_point<Clock,duration>;

  duration tt_integral = 0s;
  duration tt_free_t = 0s;
  duration tt_sing_t = 0s;
  duration tt_doub_t = 0s;

  duration tt_freeenergy = 0s;

  
  
  for (int it = 1+gp.startstep; it <= gp.steps+gp.startstep; it ++) {
    
    // compute nl(t) given phi(t), and also free energy derivatives for the different
    //   nucleation sites.

    auto current_time = Clock::now();
    
    free_energy_derivs = integrator.nonlinear(nonlinear,phi,X_is,free_energy); 

    tt_freeenergy += Clock::now() - current_time;
    
    
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

      current_time = Clock::now();




      for (int index = 0; index < beads_per_pmer.at(pmer->number); index ++) {
	transfer_free_energy(pmer->dFdX_is[index],free_energy_derivs,
			     pmer->nuc_index+index+stationary_nuc_count);

      }

      pmer->single_step(t,gp.dt,pmer->dFdX_is);

      for (int index = 0; index < beads_per_pmer.at(pmer->number); index ++) {
	transfer_nucleation_site(X_is[pmer->nuc_index+index+stationary_nuc_count],
				 pmer->atoms[pmer->nuc_beads[index]].R,
				 gp.realspace.get_Lx(),gp.realspace.get_Ly(),
				 gp.realspace.get_Lz());
      }
      

      tt_free_t += Clock::now() - current_time;      
      
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
      
      current_time = Clock::now();


      for (int index = 0; index < beads_per_pmer.at(pmer->number); index ++) {
	transfer_free_energy(pmer->dFdX_is[index],free_energy_derivs,
			     pmer->nuc_index+index+stationary_nuc_count);

      }

      pmer->single_step(t,gp.dt,pmer->dFdX_is);

      for (int index = 0; index < beads_per_pmer.at(pmer->number); index ++) {
	transfer_nucleation_site(X_is[pmer->nuc_index+index+stationary_nuc_count],
				 pmer->atoms[pmer->nuc_beads[index]].R,
				 gp.realspace.get_Lx(),gp.realspace.get_Ly(),
				 gp.realspace.get_Lz());
      }

      tt_sing_t += Clock::now() - current_time;
      
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


      current_time = Clock::now();

      for (int index = 0; index < beads_per_pmer.at(pmer->number); index ++) {
	transfer_free_energy(pmer->dFdX_is[index],free_energy_derivs,
			     pmer->nuc_index+index+stationary_nuc_count);

      }

      pmer->single_step(t,gp.dt,pmer->dFdX_is);

      for (int index = 0; index < beads_per_pmer.at(pmer->number); index ++) {
	transfer_nucleation_site(X_is[pmer->nuc_index+index+stationary_nuc_count],
				 pmer->atoms[pmer->nuc_beads[index]].R,
				 gp.realspace.get_Lx(),gp.realspace.get_Ly(),
				 gp.realspace.get_Lz());
      }

      tt_doub_t += Clock::now() - current_time;
      
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
    
    current_time = Clock::now();
    timestep.update(t,integrator);
    tt_integral = Clock::now() - current_time;
    
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
    broadcast_X_is(stationary_nuc_count,beads_per_pmer,X_is,gp.mpi_size,gp.comm);
    
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


  const double tot_integral = tt_integral/chronoitvl;
  const double tot_free_t = tt_free_t/chronoitvl;
  const double tot_sing_t = tt_sing_t/chronoitvl;
  const double tot_doub_t = tt_doub_t/chronoitvl;
  const double tot_freeenergy = tt_freeenergy/chronoitvl;

  std::cout << "total integration time on process " << gp.id << " is " << tot_integral
	    << std::endl;


  std::cout << "total free tether time on process " << gp.id << " is " << tot_free_t
	    << std::endl;

  std::cout << "total single tether time on process " << gp.id << " is " << tot_sing_t
	    << std::endl;

  std::cout << "total double time on process " << gp.id << " is " << tot_doub_t
	    << std::endl;


  std::cout << "total free energy time on process " << gp.id << " is " << tot_freeenergy
	    << std::endl;

  
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
