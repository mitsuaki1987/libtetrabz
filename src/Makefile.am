MAINTAINERCLEANFILES = Makefile.in

include_HEADERS = libtetrabz.h
nodist_include_HEADERS = libtetrabz.$(FC_MODEXT)
lib_LTLIBRARIES = libtetrabz.la

libtetrabz_la_SOURCES = libtetrabz.F90 libtetrabz_polcmplx_mod.F90 libtetrabz_fermigr_mod.F90 \
						libtetrabz_polstat_mod.F90 libtetrabz_dbldelta_mod.F90 libtetrabz_dblstep_mod.F90 \
						libtetrabz_dos_mod.F90 libtetrabz_occ_mod.F90 libtetrabz_common.F90

# Clean up modules.
clean-local: clean-modules
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)

libtetrabz.lo : libtetrabz.F90 libtetrabz_polcmplx_mod.lo libtetrabz_fermigr_mod.lo libtetrabz_polstat_mod.lo libtetrabz_dbldelta_mod.lo libtetrabz_dblstep_mod.lo libtetrabz_dos_mod.lo libtetrabz_occ_mod.lo
libtetrabz_common.lo : libtetrabz_common.F90
libtetrabz_dbldelta_mod.lo : libtetrabz_dbldelta_mod.F90 libtetrabz_common.lo
libtetrabz_dblstep_mod.lo : libtetrabz_dblstep_mod.F90 libtetrabz_common.lo
libtetrabz_dos_mod.lo : libtetrabz_dos_mod.F90 libtetrabz_common.lo
libtetrabz_fermigr_mod.lo : libtetrabz_fermigr_mod.F90 libtetrabz_common.lo
libtetrabz_occ_mod.lo : libtetrabz_occ_mod.F90 libtetrabz_common.lo
libtetrabz_polcmplx_mod.lo : libtetrabz_polcmplx_mod.F90 libtetrabz_common.lo
libtetrabz_polstat_mod.lo : libtetrabz_polstat_mod.F90 libtetrabz_common.lo
