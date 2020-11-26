cmake_minimum_required(VERSION 3.12)

project(GA_AEM VERSION 2.1
DESCRIPTION "Geoscience Australia AEM Inversion Software"
LANGUAGES CXX)

option(WITH_OMP_GARJMCMC "WITH_OMP_GARJMCMC" OFF)
option(WITH_MPI "WITH_MPI" ON)

#add object libraries for Ross's IO and data types
add_subdirectory(submodules/cpp-utils)

add_executable(
  garjmcmctdem
  src/garjmcmctdem.cpp
 )
target_link_libraries(garjmcmctdem PRIVATE file_utils general_utils)

add_executable(
  gaforwardmodeltdem
  src/gaforwardmodeltdem.cpp
)
target_link_libraries(gaforwardmodeltdem PRIVATE file_utils general_utils)

# cmake build broken for galeiallatonce because I can't find a nice way
# to search for petsc
add_executable(
  galeiallatonce
  EXCLUDE_FROM_ALL
  src/galeiallatonce.cpp
)
target_link_libraries(galeiallatonce PRIVATE file_utils general_utils)

add_executable(
  galeisbstdem
  src/galeisbstdem.cpp
)
target_link_libraries(galeisbstdem PRIVATE file_utils general_utils)

set_target_properties(
  garjmcmctdem
  gaforwardmodeltdem
  galeiallatonce
  galeisbstdem
  PROPERTIES
  CXX_STANDARD 11
  CXX_STANDARD_REQUIRED ON
)

if(WITH_OMP_GARJMCMC)
  message("Configuring with OpenMP multithreading...")
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(garjmcmctdem PUBLIC OpenMP::OpenMP_CXX)
  endif()
endif()

#configure MPI build for programs that support it
if(WITH_MPI)
  find_package(MPI)
  if(MPI_FOUND)
    message("Configuring with MPI support...")
    target_link_libraries(garjmcmctdem PRIVATE MPI::MPI_CXX)
    target_compile_definitions(garjmcmctdem PRIVATE _MPI_ENABLED)
    target_link_libraries(galeisbstdem PRIVATE MPI::MPI_CXX)
    target_compile_definitions(galeisbstdem PRIVATE _MPI_ENABLED)
    target_link_libraries(galeiallatonce PRIVATE MPI::MPI_CXX)
    target_compile_definitions(galeiallatonce PRIVATE _MPI_ENABLED)
  else()
    message("Could not find MPI on this system. Configuring without MPI support...")
  endif()
else()
  message("Build without MPI specified. Configuring without MPI support...")
endif()

target_include_directories(galeisbstdem PRIVATE submodules/eigen)
target_link_libraries(garjmcmctdem PRIVATE netcdf netcdf_c++4 fftw3)
target_link_libraries(galeisbstdem PRIVATE netcdf netcdf_c++4 fftw3)
target_link_libraries(galeiallatonce PRIVATE netcdf netcdf_c++4 fftw3)
target_link_libraries(gaforwardmodeltdem PRIVATE fftw3)

install(
  TARGETS garjmcmctdem galeisbstdem gaforwardmodeltdem
  RUNTIME DESTINATION bin
)
