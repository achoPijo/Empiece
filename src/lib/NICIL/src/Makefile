#
#--Makefile to remake the executable
#
FC      = gfortran
FFLAGS  = -fdefault-real-8  -fdefault-double-8 -fopenmp  -std=f2008 -fall-intrinsics
#FFLAGS  = -fdefault-real-8  -fdefault-double-8 -fopenmp  -std=f2008 -fall-intrinsics -g -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow,inexact,denormal -fcheck=all
#FFLAGS  = -fdefault-real-8  -fdefault-double-8 -fopenmp  -std=f2008 -fall-intrinsics -g -fbacktrace -fcheck=all
#FFLAGS = 
#FC     = ifort
#FFLAGS = -r8 -openmp
.SUFFIXES: .F90 

%.o : %.F90
	$(FC) -c $(FFLAGS) $< -o $@

eta: nicil.o  nicil_ex_eta_sup.o nicil_ex_eta.o
	$(FC) $(FFLAGS) -o nicil_ex_eta nicil.o  nicil_ex_eta_sup.o nicil_ex_eta.o
	cp nicil_ex_eta ..

sph: nicil.o  nicil_ex_sph_sup.o nicil_ex_sph.o
	$(FC) $(FFLAGS) -o nicil_ex_sph nicil.o  nicil_ex_sph_sup.o nicil_ex_sph.o
	cp nicil_ex_sph ..

clean:   
	rm -f *.o *.mod
	rm nicil_ex_eta; rm nicil_ex_sph

allclean:   
	rm -f *.o *.mod
	rm -f ../data/*.dat; 	rm -f ../data/*.log
	rm -f ../Graphs/plot.gnuplot; rm -f ../Graphs/*eps
	rm ../nicil_ex_eta; rm ../nicil_ex_sph; rm nicil_ex_eta; rm nicil_ex_sph
