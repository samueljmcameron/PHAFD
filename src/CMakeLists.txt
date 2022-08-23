add_executable(ll_w_p main.cpp run.cpp globalparams.cpp)

target_link_libraries(ll_w_p PUBLIC MPI::MPI_CXX PkgConfig::PSPDE PkgConfig::FFTW ${FFTW_MPI_LIB} PkgConfig::BEADRODPMER Eigen3::Eigen)
