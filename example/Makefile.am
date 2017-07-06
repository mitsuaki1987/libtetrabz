MAINTAINERCLEANFILES    = Makefile.in

noinst_PROGRAMS = dos.x lindhard.x dos_c.x

dos_x_SOURCES = dos.F90
lindhard_x_SOURCES = dos.F90
dos_c_x_SOURCES = dos_c.c

AM_FCFLAGS = -I$(top_builddir) -I$(top_builddir)/src
AM_CFLAGS = -I$(top_builddir) -I$(top_builddir)/src
dos_x_LDADD = -L$(top_builddir)/src -ltetrabz $(LAPACK_LIBS) $(BLAS_LIBS)
lindhard_x_LDADD = -L$(top_builddir)/src -ltetrabz $(LAPACK_LIBS) $(BLAS_LIBS)
dos_c_x_LDADD = -L$(top_builddir)/src -ltetrabz $(LAPACK_LIBS) $(BLAS_LIBS) -lm
