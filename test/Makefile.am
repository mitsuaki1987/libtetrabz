MAINTAINERCLEANFILES = Makefile.in
EXTRA_DIST = $(TESTS)

TESTS = test2_8_8.sh test2_16_8.sh test2_16_16.sh

check_PROGRAMS = test.x
test_x_SOURCES = test.F90
AM_FCFLAGS = $(FC_MODINC)$(top_builddir) $(FC_MODINC)$(top_builddir)/src 
test_x_LDADD = -L$(top_builddir)/src -ltetrabz

distclean-local:
	-rm -rf test.in test.out test.out.ref

# Clean up modules.
clean-local: clean-modules
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
