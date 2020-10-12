#####################################################################
# Command Line Options (see man page)
#  -B   make all parts from scratch   
#
# Debug Options 
#
#  -g 
#  -check all
#
# should be removed for speed tests
#
#

FC      = gfortran -g -fbacktrace -fcheck=all -Wall -Wno-tabs
FFLAGS  = -O3
CFLAGS  = -O -DLINUX
LINKER  = gfortran -g
LIBS    =
PROG    = ./scf

#####################################################################
# manually add one object file for each source file (.f or .f90)
OBJS= main.o int.o slater.o prmat.o nuclearrepulsion.o atomtype.o matrixpacker.o matrixunpacker.o echo.o setupbasis.o scfcycle.o dspev.o aufpunkt.o matrices.o mullikanpop.o electrondensity.o surfacescan.o gradient.o
#####################################################################
# compile a .o from each .f
%.o: %.f
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) -c $< -o $@

# compile a .o from each .f90
%.o: %.f90
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) -c $< -o $@
# link
$(PROG): $(OBJS)
	$(LINKER) $(OBJS) $(LIBS) -o $(PROG)

clean:
	rm -f *.o $(PROG)

