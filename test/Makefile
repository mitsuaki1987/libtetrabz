include ../make.sys

all:pdos.x tetlambda.x pdos_mpi.x tetlambda_mpi.x chi0.x occstep.x

clean:
	rm -f *.x *.a *.o *.mod

pdos.x:pdos.f90 ../src/libtetrabz.a
	$(F90) $< -L $(TOPDIR)/src/ -I $(TOPDIR)/src/ -ltetrabz $(FFLAGS) -o $@ 
	rm -f *.o *.mod

pdos_mpi.x:pdos_mpi.f90 ../src/libtetrabz_mpi.a
	$(MPIF90) $< -L $(TOPDIR)/src/ -I $(TOPDIR)/src/ -ltetrabz_mpi $(FFLAGS) -o $@ 
	rm -f *.o *.mod

tetlambda.x:tetlambda.f90 ../src/libtetrabz.a
	$(F90) $< -L $(TOPDIR)/src/ -I $(TOPDIR)/src/ -ltetrabz $(FFLAGS) -o $@ 
	rm -f *.o *.mod

tetlambda_mpi.x:tetlambda_mpi.f90 ../src/libtetrabz_mpi.a
	$(MPIF90) $< -L $(TOPDIR)/src/ -I $(TOPDIR)/src/ -ltetrabz_mpi $(FFLAGS) -o $@ 
	rm -f *.o *.mod

chi0.x:chi0.f90 ../src/libtetrabz_mpi.a
	$(MPIF90) $< -L $(TOPDIR)/src/ -I $(TOPDIR)/src/ -ltetrabz_mpi $(FFLAGS) -o $@ 
	rm -f *.o *.mod

occstep.x:occstep.f90 ../src/libtetrabz_mpi.a
	$(MPIF90) $< -L $(TOPDIR)/src/ -I $(TOPDIR)/src/ -ltetrabz_mpi $(FFLAGS) -o $@ 
	rm -f *.o *.mod

