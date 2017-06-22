include ../make.sys

all:dos.x lindhard.x dos_c.x

clean:
	rm -f *.x *.a *.o *.mod

dos.x:dos.F90 ../src/libtetrabz.a
	$(F90) $< -L ../src/ -I ../src/ -ltetrabz $(FFLAGS) -o $@ 
	rm -f *.o *.mod

lindhard.x:lindhard.F90 ../src/libtetrabz.a
	$(F90) $< -L ../src/ -I ../src/ -ltetrabz $(FFLAGS) -o $@ 
	rm -f *.o *.mod

dos_c.x:dos.c ../src/libtetrabz.a
	$(CC) $< -L ../src/ -I ../src/ -ltetrabz $(CFLAGS) -o $@ 
	rm -f *.o

