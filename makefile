# Compiler
COMPILE=gfortran -g -O3 -fcheck=all -Wall# -O3 #-g -O0 #-O3 #-g #-O3
#COMPILE=mpif90

main_2D_reconstruction : main_2D_reconstruction.o libBezier.o TimeScheme.o AdvectionProblem.o interpol.o variables.o midline_2D_modules.o
	$(COMPILE) -o  main_2D_reconstruction main_2D_reconstruction.o libBezier.o AdvectionProblem.o TimeScheme.o interpol.o variables.o midline_2D_modules.o $(LIBS)

main_3D_reconstruction : main_3D_reconstruction.o libBezier.o TimeScheme.o AdvectionProblem.o interpol.o variables.o midline_2D_modules.o
	$(COMPILE) -o  main_3D_reconstruction main_3D_reconstruction.o libBezier.o AdvectionProblem.o TimeScheme.o interpol.o variables.o midline_2D_modules.o $(LIBS)

main_2D_reconstruction.o : main_2D_reconstruction.f90 libBezier.mod TimeScheme.mod AdvectionProblem.mod interpol.mod variables.mod midline_2D_modules.mod
	$(COMPILE) -c  main_2D_reconstruction.f90

main_3D_reconstruction.o : main_3D_reconstruction.f90 libBezier.mod TimeScheme.mod AdvectionProblem.mod interpol.mod variables.mod midline_2D_modules.mod
	$(COMPILE) -c  main_3D_reconstruction.f90

AdvectionProblem.o AdvectionProblem.mod : AdvectionProblem.f90 interpol.mod variables.mod doubler.mod 
	$(COMPILE) -c  AdvectionProblem.f90

TimeScheme.o TimeScheme.mod : TimeScheme.f90 AdvectionProblem.mod  
	$(COMPILE) -c  TimeScheme.f90

interpol.o interpol.mod: interpol.f90 doubler.mod
	        $(COMPILE) -c   interpol.f90

libBezier.o libBezier.mod: libBezier.f90 AdvectionProblem.mod doubler.mod
	        $(COMPILE) -c   libBezier.f90

midline_2D_modules.o midline_2D_modules.mod: midline_2D_modules.f90 AdvectionProblem.mod doubler.mod
	        $(COMPILE) -c   midline_2D_modules.f90

variables.o variables.mod: variables.f90 doubler.mod
	$(COMPILE) -c  variables.f90 

doubler.o doubler.mod: doubler.f90
	$(COMPILE) -c   doubler.f90 

# Remove object files
clean :
	rm -f *.o *.mod u *.vtk deformation3D_para deformation3D deformation2D advectionMain

