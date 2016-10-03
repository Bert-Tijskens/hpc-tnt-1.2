# ============================================================================
# Name        : Makefile
# Author      : 
# Version     :
# Copyright   : Your copyright notice
# Description : Makefile for Hello World in Fortran
# ============================================================================

.PHONY: all clean

# ${fortran_compiler} is defined in project properties/Fortran Build/environment
FORTRAN_COMPILER = ${fortran_compiler}

ADVISOR_OPTIONS = -vec -simd -g -qopt-report=5  -xHost

all: src/ppmd01.f90
	$(FORTRAN_COMPILER) -O2 -g $(ADVISOR_OPTIONS) \
	    -o ppmd01 \
	    src/prec.f90 \
		src/wprec.f90 \
		src/md.f90 \
		src/util_random.f90 \
		src/ppmd01.f90

clean:
	rm -f bin/ppmd01 *.mod

#	ifort -O2 -g -vec -simd -g -qopt-report=5  -xHost -S src/prec.f90 src/wprec.f90 src/md.f90 src/util_random.f90 src/ppmd01.f90