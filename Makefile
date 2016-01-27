
all:libs examples

libs:
	cd src; make

examples:
	cd example; make

clean:
	cd src; make clean
	cd example; make clean

