# ==============================================================================
# makefile for DAISY code
# Version 1.0.0
# 2015-04-08
# ==============================================================================

program = DAISY.linux
source = $(wildcard *.f90)
objects = $(source:.f90=.o)

prefix = ../lib-linux
prehdf5 = /apps/compile/hdf5/1.8.16-intel-16.0.0
path_inc = 
path_lib = 

# include
inc_hdf5    = -I$(prehdf5)/include 

inc_collect = -I$(prefix)/collection/include 
inc_xml     = -I$(prefix)/fox/include 
inc_lapack  = -I$(prefix)/lapack/include 
inc_refprop = -I$(prefix)/refprop/include 
path_inc += $(inc_collect)
path_inc += $(inc_xml)
path_inc += $(inc_hdf5)
path_inc += $(inc_lapack)
path_inc += $(inc_refprop)

# library
lib_hdf5    = -L$(prehdf5)/lib       -lhdf5 -lhdf5_fortran -lhdf5_hl -lhdf5hl_fortran

lib_collect = -L$(prefix)/collection/lib -lcollection
lib_xml     = -L$(prefix)/fox/lib        -lfox_common -lfox_dom -lfox_fsys -lfox_sax -lfox_utils -lfox_wxml
lib_lapack  = -L$(prefix)/lapack/lib     -ltmglib -llapack -lrefblas
lib_refprop = -L$(prefix)/refprop/lib    -lrefprop
path_lib += $(lib_collect)
path_lib += $(lib_xml)
path_lib += $(lib_hdf5)
path_lib += $(lib_lapack)
path_lib += $(lib_refprop)

# compiler
COMPILER = intel
ifeq ($(COMPILER),intel)
  F90 = ifort
  F90FLAGS = -O2 -openmp
  LDFLAGS = -openmp
endif

# ==============================================================================
pre:
	@cp driver/*.f90 .
	@cp frame/*.f90 .
	@cp header/*.f90 .
	@cp inout/*.f90 .
	@cp NK2solver/*.f90 .
	@cp solver/*.f90 .
	@cp thermal/*.f90 .
	@dos2unix *.f90    

build: $(program)
$(program): $(objects)
	$(F90) $(objects) $(path_lib) $(LDFLAGS) -o $@
    
clean: 
	@rm -f *.o *.mod $(program)
    
post: 
	@rm -f *.o *.mod *.f90
    
# ==============================================================================
.SUFFIXES: .f90 .o
.PHONY: pre build clean post

%.o: %.f90
	$(F90) $(F90FLAGS) $(path_inc) -c $< 
    
include dependence.inc
