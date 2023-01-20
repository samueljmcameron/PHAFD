add_executable(ll_w_p main.cpp atom.cpp domain.cpp read_atoms.cpp utility.cpp input.cpp group.cpp fixatom.cpp fixatom_semiflexible.cpp comm_brick.cpp nbin.cpp nbin_standard.cpp neigh_list.cpp neighbor.cpp npair.cpp npair_grid_bin.cpp npair_half_bin_newton.cpp nstencil.cpp nstencil_full_bin_3d.cpp nstencil_half_bin_3d.cpp my_page.cpp pair.cpp pair_gridatom_gaussian.cpp pair_lj_cut.cpp)

target_link_libraries(ll_w_p PUBLIC MPI::MPI_CXX PkgConfig::PSPDE PkgConfig::FFTW ${FFTW_MPI_LIB} PkgConfig::BEADRODPMER Eigen3::Eigen)
