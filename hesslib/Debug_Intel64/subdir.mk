################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../num_types.o \
../numerical_utilities.o \
../symmetry_module.o \
../vector_matrix_utilities.o 

F90_SRCS += \
../compare_structures.f90 \
../num_types.f90 \
../numerical_utilities.f90 \
../symmetry_module.f90 \
../vector_matrix_utilities.f90 

OBJS += \
./compare_structures.o \
./num_types.o \
./numerical_utilities.o \
./symmetry_module.o \
./vector_matrix_utilities.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: Intel(R) Intel(R) 64 Fortran Compiler'
	ifort -g -O0 -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


