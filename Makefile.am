ACLOCAL_AMFLAGS = -I config/m4
SUBDIRS = src test example
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 configure
EXTRA_DIST = docs

doc:
	cd docs/ja/;make latexpdf; make html
	cd docs/en/;make latexpdf; make html

pngdoc:
	cd docs/ja/;sed -i -e "s/mathjax/pngmath/g" conf.py; make latexpdf; make html
	cd docs/en/;sed -i -e "s/mathjax/pngmath/g" conf.py; make latexpdf; make html

doc-clean:
	cd docs/ja/; make clean; sed -i -e "s/pngmath/mathjax/g" conf.py
	cd docs/en/; make clean; sed -i -e "s/pngmath/mathjax/g" conf.py
