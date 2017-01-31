# Set up conditional compilation so that one Makefile can handle
# multiple platforms/compilers. Branch according to an environmental
# variable F90. I wish someone would show me a better way of doing this.

ifeq (${F90},gfortran)  # free fortran compiler
  FFLAGS =  -g -fbounds-check -fPIC -pedantic -Wall -ffree-line-length-none
  FOUND = true
endif

ifeq (${F90},ifort)  # Intel compiler
  ifeq (${DEBUG},false)
     FFLAGS =  -O3 -fPIC
     FOUND = true
  endif
  ifeq (${DEBUG},true)
     FFLAGS =  -g -error-limit 7 -traceback -check bounds -check uninit -check pointers -fp-stack-check  -warn -heap-arrays -fpic
     FOUND = true
  endif
  ifeq (${DEBUG},)
     FFLAGS =  -g -error-limit 7 -traceback -check bounds -warn -heap-arrays -fpic
     FOUND = true
  endif

endif

ifeq (${F90},ifc)  # Intel compiler
  F90 =  /opt/intel/fc/9.1.040/bin/ifort
  FFLAGS =  -g -fPIC -error-limit 7 -traceback -check bounds -warn -heap-arrays
  FOUND = true
endif
ifeq (${F90},xlf90) # IBM compiler
  FFLAGS = -g -C -fPIC -qsuffix=f=f90  
  FOUND = true
endif
ifeq (${F90},f95) # Absoft PPC compiler
  FFLAGS = -Rb -Rc  -z2 -ea -fPIC
# B80 entry into each subroutine  ; Rb bounds; Rc array conformance
  FOUND = true
endif
ifeq (${F90},)  # If the variable isn't set, make sure compilation bombs
error:
	echo Error environmental variable "F90" not set!
endif
ifneq (${FOUND},true) # If it's set but no matching flags then bomb
error:	
	echo Error: makefile doesn\'t have flags for this compiler
endif

SRC =  num_types.f90 numerical_utilities.f90 vector_matrix_utilities.f90 symmetry_module.f90
OBJS = ${SRC:.f90=.o}


.SUFFIXES :  
.SUFFIXES :  .f .f90 .f95 .o

all: hesslib.a hesslib.so

hesslib.a: num_types.o numerical_utilities.o vector_matrix_utilities.o symmetry_module.o compare_structures.o #wrapperlib.o
	ar ru $@ $?
	ranlib  $@

hesslib.so:  num_types.o numerical_utilities.o vector_matrix_utilities.o  symmetry_module.o compare_structures.o #wrapperlib.o
	ifort -shared -fpic -o hesslib.so num_types.f90  numerical_utilities.f90 vector_matrix_utilities.f90 \
	          symmetry_module.f90 compare_structures.f90
#	          symmetry_module.f90 compare_structures.f90 wrapperlib.f90

#vec_mat_utilities.o: 
	

num_types.o : num_types.f90
	${F90} ${FFLAGS} -c $<

numerical_utilities.o : numerical_utilities.f90
	${F90} ${FFLAGS} -c $<

vector_matrix_utilities.o : vector_matrix_utilities.f90
	${F90} ${FFLAGS} -c $<

symmetry_module.o : ../hesslib/symmetry_module.f90 
	${F90} ${FFLAGS} -c $<

compare_structures.o : ../hesslib/compare_structures.f90
	${F90} ${FFLAGS} -c $<	
	
#wrapperlib.o : wrapperlib.f90
#	${F90} ${FFLAGS} -c $<

CLEAN  = *.o *.mod *.a
clean : 
	rm -f ${CLEAN}

remake:
	make clean
	make

clobber: 
	rm -f  ${CLEAN}  *~ \#*
	make
