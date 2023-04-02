cmake_minimum_required(VERSION 3.14)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


project(phafd VERSION 1.0 DESCRIPTION "Parallel Hybrid-Atom Field Dynamics"
        LANGUAGES CXX)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")


include(CheckIPOSupported)
check_ipo_supported(RESULT supported OUTPUT error)

if( supported )
    message(STATUS "IPO / LTO enabled")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
	message(STATUS "IPO / LTO not supported: <${error}>")
endif()


find_package(MPI REQUIRED)
add_definitions(-DOMPI_SKIP_MPICXX)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

find_package(PkgConfig REQUIRED)

set(PKG_CHECK_MODULES_STATIC)

set(CMAKE_FIND_LIBRARY_SUFFIXES .a)

pkg_check_modules(FFTW-MPI IMPORTED_TARGET QUIET fftw3mpi)
if (NOT FFTW-MPI_FOUND)
  message(FATAL_ERROR "Need to specify where to find fftw_arr in CMAKE_PREFIX_PATH")
endif()

pkg_check_modules(FFTWARR IMPORTED_TARGET REQUIRED fftw_arr)

pkg_check_modules(BEADRODPMER IMPORTED_TARGET REQUIRED beadrodpmer)


add_subdirectory(src)



MESSAGE(STATUS "Compiler flags:" ${CMAKE_CXX_COMPILE_FLAGS})
MESSAGE(STATUS "Compiler cxx debug flags:" ${CMAKE_CXX_FLAGS_DEBUG})
MESSAGE(STATUS "Compiler cxx release flags:" ${CMAKE_CXX_FLAGS_RELEASE})
MESSAGE(STATUS "Compiler cxx min size flags:" ${CMAKE_CXX_FLAGS_MINSIZEREL})
MESSAGE(STATUS "Compiler cxx flags:" ${CMAKE_CXX_FLAGS})
