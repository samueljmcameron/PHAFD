cmake_minimum_required(VERSION 3.14)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


project(liqliq_with_polymers VERSION 1.0 DESCRIPTION "Model B Dynamics"
        LANGUAGES CXX)

set(ENV{PKG_CONFIG_PATH} "${PKG_FFTW_MPI_LIB_DIR}/pkgconfig")

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

# this is from alfC answer at url https://stackoverflow.com/questions/44277114/how-to-use-fftw-library-in-cmake
find_package(PkgConfig REQUIRED)

set(PKG_CHECK_MODULES_STATIC)

set(CMAKE_FIND_LIBRARY_SUFFIXES .a)

pkg_check_modules(FFTW IMPORTED_TARGET REQUIRED fftw3)

if( NOT FFTW_ROOT AND DEFINED ENV{FFTWDIR} )
    set( FFTW_ROOT $ENV{FFTWDIR} )
endif()

find_library(FFTW_MPI_LIB NAMES "libfftw3_mpi.a" PATHS ${PKG_FFTW_MPI_LIB_DIR})

if( FFTW_MPI_LIB)
    set(FFTW_MPI_LIB_FOUND TRUE)
    set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_MPI_LIB})
    add_library(FFTW::MPI INTERFACE IMPORTED)
    set_target_properties(FFTW::MPI PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}" INTERFACE_LINK_LIBRARIES "${FFTW_MPI_LIB}")
else()
    set(FFTW_MPI_LIB_FOUND FALSE)
endif()

pkg_check_modules(PSPDE IMPORTED_TARGET REQUIRED ps_pde)

pkg_check_modules(BEADRODPMER IMPORTED_TARGET REQUIRED beadrodpmer)


add_subdirectory(src)



MESSAGE(STATUS "Compiler flags:" ${CMAKE_CXX_COMPILE_FLAGS})
MESSAGE(STATUS "Compiler cxx debug flags:" ${CMAKE_CXX_FLAGS_DEBUG})
MESSAGE(STATUS "Compiler cxx release flags:" ${CMAKE_CXX_FLAGS_RELEASE})
MESSAGE(STATUS "Compiler cxx min size flags:" ${CMAKE_CXX_FLAGS_MINSIZEREL})
MESSAGE(STATUS "Compiler cxx flags:" ${CMAKE_CXX_FLAGS})
