include ../make.sys

all:pdos.x tetlambda.x

clean:
	rm -f *.x *.a *.o *.mod

pdos.x:pdos.f90 ../src/libtetrabz.a
	$(F90) $< -L $(TOPDIR)/src/ -I $(TOPDIR)/src/ -ltetrabz $(FFLAGS) -o $@ 

tetlambda.x:tetlambda.f90 ../src/libtetrabz.a
	$(F90) $< -L $(TOPDIR)/src/ -I $(TOPDIR)/src/ -ltetrabz $(FFLAGS) -o $@ 

