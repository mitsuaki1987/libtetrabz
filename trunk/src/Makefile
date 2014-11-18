include ../make.sys

all:libtetrabz.a libtetrabz_mpi.a

clean:
	rm -f *.a *.o *.mod

libtetrabz.a:libtetra_common.o libtetra.o
	ar cr $@ libtetra_common.o libtetra.o

libtetrabz_mpi.a:libtetra_common.o libtetra_mpi.o
	ar cr $@ libtetra_common.o libtetra_mpi.o

libtetra.o:libtetra_common.o libtetra.f90
	$(F90) libtetra.f90 $(FFLAGS) -c

libtetra_mpi.o:libtetra_common.o libtetra_mpi.f90
	$(MPIF90) libtetra_mpi.f90 $(FFLAGS) -c

libtetra_common.o:libtetra_common.f90
	$(F90) libtetra_common.f90 $(FFLAGS) -c
