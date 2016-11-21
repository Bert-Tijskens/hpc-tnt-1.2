# ============================================================================
# Name        : Makefile
# Author      : 
# Version     :
# Copyright   : Your copyright notice
# Description : Makefile for pyMDFortran 
# ============================================================================

.PHONY: all clean

# Change this line if you are using a different Fortran compiler
FORTRAN_COMPILER = ifort

all: src/testMDFortran.f90 src/md.f90
	$(FORTRAN_COMPILER) -O3 -xHost -vec -simd -g  -qopt-report=5 -fsource-asm\
		-o bin/testMDFortran \
		src/md.f90			 \
		src/testMDFortran.f90
	f2py -c src/md.f90 -m pyMDFortran --f90flags="-g" --arch="-xHost" --opt="-O3 -qopt-report=5 -vec -simd" --fcompiler=intelem
	
clean:
	rm -f bin/testMDFortran *.mod *.so
