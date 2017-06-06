include ../make.sys

all:libtetrabz.a libtetrabz_mpi.a

clean:
	rm -f *.a *.o *.mod

libtetrabz.a:libtetrabz_common.o libtetrabz.o
	ar cr $@ libtetrabz_common.o libtetrabz.o

libtetrabz_mpi.a:libtetrabz_common.o libtetrabz_mpi.o
	ar cr $@ libtetrabz_common.o libtetrabz_mpi.o

libtetrabz.o:libtetrabz_common.o libtetrabz.F90
	$(F90) libtetrabz.F90 $(FFLAGS) -c

libtetrabz_mpi.o:libtetrabz_common.o libtetrabz_mpi.F90
	$(MPIF90) libtetrabz_mpi.F90 $(FFLAGS) -c

libtetrabz_common.o:libtetrabz_common.F90
	$(F90) libtetrabz_common.F90 $(FFLAGS) -c
