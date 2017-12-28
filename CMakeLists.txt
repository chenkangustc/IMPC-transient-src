# ==============================================================================
# cmake file for DAISY code
# Version 1.0.0
# 2015-04-08
# ==============================================================================

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(DAISY Fortran)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)

set(prefix  "/home/mthe/software/DAISY_ifort/external")
set(path_inc  "")
set(path_lib  "")

# include
set(inc_collect "-I${prefix}/collection/include")
set(inc_xml     "-I${prefix}/FoX/include") 
set(inc_hdf5    "-I${prefix}/HDF5/include") 
set(inc_lapack  "-I${prefix}/Lapack/include") 
set(inc_refprop "-I${prefix}/REFPROP/include") 
set(path_inc "${inc_collect}  ${path_inc}")
set(path_inc "${inc_xml}      ${path_inc}")
set(path_inc "${inc_hdf5}     ${path_inc}")
set(path_inc "${inc_lapack}   ${path_inc}")
set(path_inc "${inc_refprop}  ${path_inc}")

# library
set(lib_collect "-I${prefix}/collection/lib -lcollection ")
set(lib_xml     "-I${prefix}/FoX/lib        -lxml") 
set(lib_hdf5    "-I${prefix}/HDF5/lib       -lhdf5 -lhdf5_fortran -lhdf5_hl -lhdf5hl_fortran") 
set(lib_lapack  "-I${prefix}/Lapack/lib     -ltmglib -llapack -lrefblas") 
set(lib_refprop "-I${prefix}/REFPROP/lib    -lrefprop") 
set(path_lib "${lib_collect}  ${path_lib}")
set(path_lib "${lib_xml}      ${path_lib}")
set(path_lib "${lib_hdf5}     ${path_lib}")
set(path_lib "${lib_lapack}   ${path_lib}")
set(path_lib "${lib_refprop}  ${path_lib}")

set(extra_lib  "")
set(extra_lib  ${extra_lib} )

# compiler
if(CMAKE_Fortran_COMPILER_ID STREQUAL "intel")
    set(f90flags  "-o2")
    set(ldflags  "-openmp")
endif()

# ==============================================================================
set(program "DAISY")
file(GLOB source driver/*.f90 frame/*.f90 header/*.f90 inout/*.f90 NK2solver/*.f90 solver/*.f90 thermal/*.f90)

add_executable(${program} ${source})
target_link_libraries(${program} ${libraries})
set_target_properties(${program} PROPERTIES
  COMPILE_FLAGS "${f90flags}"
  LINK_FLAGS "${ldflags}")
