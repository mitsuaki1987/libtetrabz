include ../make.sys

all:test.x

clean:
	rm -f *.x *.a *.o *.mod

test.x:test.F90 ../src/libtetrabz.a
	$(F90) $< -L ../src/ -I ../src/ -ltetrabz $(FFLAGS) -o $@ 


