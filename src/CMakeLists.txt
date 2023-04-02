add_executable(phafd main.cpp phafd.cpp atom.cpp domain.cpp grid.cpp comm_brick.cpp neighbor.cpp group.cpp input.cpp nbin.cpp nbin_standard.cpp neigh_list.cpp  npair.cpp npair_grid_bin.cpp npair_half_bin_newton.cpp nstencil.cpp nstencil_full_bin_3d.cpp nstencil_half_bin_3d.cpp my_page.cpp utility.cpp read_atoms.cpp fix.cpp fixatom_semiflexible.cpp fixgrid_conjugate.cpp pair.cpp pair_gridatom_gaussian.cpp pair_gridatom_lj_ish.cpp pair_lj_cut.cpp pair_harmonic_cut.cpp fixgrid_floryhuggins.cpp compute.cpp compute_complex.cpp dump.cpp fixgrid_ave.cpp integrate.cpp compute_pair.cpp fixgrid_gradphi.cpp fixatom_drag.cpp pair_gridatom_lj_ish_linear.cpp read_dump.cpp conjugate.cpp conjugate_volfrac.cpp smatrix.cpp)

target_link_libraries(phafd PUBLIC MPI::MPI_CXX PkgConfig::FFTWARR PkgConfig::FFTW-MPI PkgConfig::BEADRODPMER Eigen3::Eigen)

target_include_directories(phafd PUBLIC PkgConfig::FFTW-MPI)
