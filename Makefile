
all:
	cd src; make
	cd example; make
	cd test; make

clean:
	cd test; make clean
	cd example; make clean
	cd src; make clean

