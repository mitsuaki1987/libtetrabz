MAINTAINERCLEANFILES = Makefile.in
noinst_PROGRAMS = test.x

test_x_SOURCES = test.F90

AM_FCFLAGS = -I$(top_builddir) -I$(top_builddir)/src 
test_x_LDADD = -L$(top_builddir)/src -ltetrabz
