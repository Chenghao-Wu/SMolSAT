cmake_minimum_required(VERSION 3.10)

# set the project name
project(SMolSAT)

# set install directory
set(INSTALL_PY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/SMolSAT)

# find python environment
find_package(Python COMPONENTS Interpreter Development REQUIRED)

# add subdirectory of pybind11 and qvectors
add_subdirectory(./extern/pybind11)

# set qvector directory for qvectors
add_subdirectory(./extern/qvectors)
set(WAVEVECTORS3D \"${CMAKE_CURRENT_SOURCE_DIR}/extern/qvectors/qvectors3d/qvector\")
set(WAVEVECTORS2D \"${CMAKE_CURRENT_SOURCE_DIR}/extern/qvectors/qvectors2d/qvector\")
set(WAVEVECTORS1D \"${CMAKE_CURRENT_SOURCE_DIR}/extern/qvectors/qvectors1d/qvector\")
add_definitions(-DWV3D=${WAVEVECTORS3D})
add_definitions(-DWV2D=${WAVEVECTORS2D})
add_definitions(-DWV1D=${WAVEVECTORS1D})

# set EIGEN directory
set(EIGEN_DIR "./extern/eigen-3.3.9")
include_directories(${EIGEN_DIR})

# set source files
set(_pmdat_sources ./SMolSAT/coordinate.cpp
                    ./SMolSAT/trajectory.cpp
                    ./SMolSAT/atom_trajectory.cpp
                    ./SMolSAT/molecule.cpp
                    ./SMolSAT/error.cpp
                    ./SMolSAT/multibody_set.cpp
                    ./SMolSAT/progress.cpp
                    ./SMolSAT/system.cpp
                    ./SMolSAT/multibody.cpp
                    ./SMolSAT/trajmath.cpp
                    ./SMolSAT/tokenize.cpp
                    ./SMolSAT/trajectory_set.cpp
                    ./SMolSAT/multibody_analysis.cpp
                    ./SMolSAT/multibody_list.cpp
                    ./SMolSAT/trajectory_list.cpp
                    ./SMolSAT/static_trajectory_list.cpp
                    ./SMolSAT/analysis_base.cpp
                    ./SMolSAT/analysis_onetime.cpp
                    ./SMolSAT/correlation_2d.cpp
                    ./SMolSAT/boolean_list.cpp
                    ./SMolSAT/value_list.cpp
                    ./SMolSAT/trajectories.cpp
                    ./SMolSAT/mean_square_displacement.cpp
                    ./SMolSAT/displacement_list.cpp
                    ./SMolSAT/bond_autocorrelation_function.cpp
                    ./SMolSAT/rgtensor.cpp
                    ./SMolSAT/rgtensor_stats.cpp
                    ./SMolSAT/wave_vectors.cpp
                    ./SMolSAT/structure_factor.cpp
                    ./SMolSAT/radial_distribution_function.cpp
                    ./SMolSAT/non_gaussian_parameter.cpp
                    ./SMolSAT/incoherent_scattering_function.cpp
                    ./SMolSAT/space-time_correlation_function.cpp
                    ./SMolSAT/trajectory_list_bins.cpp
                    ./SMolSAT/van_hove_self.cpp
                    ./SMolSAT/van_hove_distinct.cpp
                    ./SMolSAT/pmdat.cpp)

# set header files
set(_pmdat_headers ./SMolSAT/coordinate.h
                    ./SMolSAT/trajectory.h
                    ./SMolSAT/atom_trajectory.h
                    ./SMolSAT/molecule.h
                    ./SMolSAT/error.h
                    ./SMolSAT/multibody_set.h
                    ./SMolSAT/progress.h
                    ./SMolSAT/system.h
                    ./SMolSAT/multibody.h
                    ./SMolSAT/trajmath.h
                    ./SMolSAT/vector_map.h
                    ./SMolSAT/tokenize.h
                    ./SMolSAT/trajectory_set.h
                    ./SMolSAT/trajectory_list.h
                    ./SMolSAT/static_trajectory_list.h
                    ./SMolSAT/multibody_analysis.h
                    ./SMolSAT/multibody_list.h
                    ./SMolSAT/analysis_base.h
                    ./SMolSAT/analysis_onetime.h
                    ./SMolSAT/space-time_correlation_function.h
                    ./SMolSAT/trajectory_list_bins.h
                    ./SMolSAT/correlation_2d.h
                    ./SMolSAT/boolean_list.h
                    ./SMolSAT/value_list.h
                    ./SMolSAT/trajectories.h
                    ./SMolSAT/mean_square_displacement.h
                    ./SMolSAT/displacement_list.h
                    ./SMolSAT/bond_autocorrelation_function.h
                    ./SMolSAT/rgtensor.h
                    ./SMolSAT/rgtensor_stats.h
                    ./SMolSAT/wave_vectors.h
                    ./SMolSAT/structure_factor.h
                    ./SMolSAT/radial_distribution_function.h
                    ./SMolSAT/non_gaussian_parameter.h
                    ./SMolSAT/incoherent_scattering_function.h
                    ./SMolSAT/van_hove_self.h
                    ./SMolSAT/van_hove_distinct.h
                    ./SMolSAT/xdrfile.h
                    ./SMolSAT/xdrfile_xtc.h
                    )

# add pybind11 module
pybind11_add_module(_SMolSAT ./SMolSAT/module.cpp ${_pmdat_sources})

add_executable(SMolSAT ${_pmdat_sources} ${_pmdat_headers})

target_link_libraries(SMolSAT PRIVATE pybind11::embed)

install(
  TARGETS
    _SMolSAT
  LIBRARY
      DESTINATION ${INSTALL_PY_DIR}
)