
all:libs examples

libs:
	cd src; make

examples:
	cd example; make

tests:
	cd test; make

clean:
	cd src; make clean
	cd example; make clean
	cd test; make clean
