MAINTAINERCLEANFILES = Makefile.in

include_HEADERS = libtetrabz.h libtetrabz.mod
lib_LTLIBRARIES = libtetrabz.la

libtetrabz_la_SOURCES = libtetrabz.F90 libtetrabz_polcmplx_mod.F90 libtetrabz_fermigr_mod.F90 \
						libtetrabz_polstat_mod.F90 libtetrabz_dbldelta_mod.F90 libtetrabz_dblstep_mod.F90 \
						libtetrabz_dos_mod.F90 libtetrabz_occ_mod.F90 libtetrabz_common.F90

libtetrabz.lo : libtetrabz.F90 libtetrabz_polcmplx_mod.lo libtetrabz_fermigr_mod.lo libtetrabz_polstat_mod.lo libtetrabz_dbldelta_mod.lo libtetrabz_dblstep_mod.lo libtetrabz_dos_mod.lo libtetrabz_occ_mod.lo
libtetrabz_common.lo : libtetrabz_common.F90
libtetrabz_dbldelta_mod.lo : libtetrabz_dbldelta_mod.F90 libtetrabz_common.lo
libtetrabz_dblstep_mod.lo : libtetrabz_dblstep_mod.F90 libtetrabz_common.lo
libtetrabz_dos_mod.lo : libtetrabz_dos_mod.F90 libtetrabz_common.lo
libtetrabz_fermigr_mod.lo : libtetrabz_fermigr_mod.F90 libtetrabz_common.lo
libtetrabz_occ_mod.lo : libtetrabz_occ_mod.F90 libtetrabz_common.lo
libtetrabz_polcmplx_mod.lo : libtetrabz_polcmplx_mod.F90 libtetrabz_common.lo
libtetrabz_polstat_mod.lo : libtetrabz_polstat_mod.F90 libtetrabz_common.lo
