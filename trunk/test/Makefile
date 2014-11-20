include ../make.sys

all:pdos.x tetlambda.x pdos_mpi.x tetlambda_mpi.x 

clean:
	rm -f *.x *.a *.o *.mod

pdos.x:pdos.f90 ../src/libtetrabz.a
	$(F90) $< -L $(TOPDIR)/src/ -I $(TOPDIR)/src/ -ltetrabz $(FFLAGS) -o $@ 

pdos_mpi.x:pdos_mpi.f90 ../src/libtetrabz_mpi.a
	$(MPIF90) $< -L $(TOPDIR)/src/ -I $(TOPDIR)/src/ -ltetrabz_mpi $(FFLAGS) -o $@ 

tetlambda.x:tetlambda.f90 ../src/libtetrabz.a
	$(F90) $< -L $(TOPDIR)/src/ -I $(TOPDIR)/src/ -ltetrabz $(FFLAGS) -o $@ 

tetlambda_mpi.x:tetlambda_mpi.f90 ../src/libtetrabz_mpi.a
	$(MPIF90) $< -L $(TOPDIR)/src/ -I $(TOPDIR)/src/ -ltetrabz_mpi $(FFLAGS) -o $@ 

