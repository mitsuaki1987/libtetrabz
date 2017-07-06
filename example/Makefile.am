MAINTAINERCLEANFILES    = Makefile.in

noinst_PROGRAMS = dos.x lindhard.x dos_c.x

dos_x_SOURCES = dos.F90
lindhard_x_SOURCES = lindhard.F90
dos_c_x_SOURCES = dos_c.c

AM_FCFLAGS = $(FC_MODINC)$(top_builddir) $(FC_MODINC)$(top_builddir)/src
AM_CFLAGS = -I$(top_builddir) -I$(top_srcdir)/src -I$(top_builddir)/src
dos_x_LDADD = -L$(top_builddir)/src -ltetrabz $(LAPACK_LIBS) $(BLAS_LIBS)
lindhard_x_LDADD = -L$(top_builddir)/src -ltetrabz $(LAPACK_LIBS) $(BLAS_LIBS)
dos_c_x_LDADD = -L$(top_builddir)/src -ltetrabz $(LAPACK_LIBS) $(BLAS_LIBS) -lm

# Clean up modules.
clean-local: clean-modules
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
