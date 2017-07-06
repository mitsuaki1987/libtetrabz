MAINTAINERCLEANFILES = Makefile.in

include_HEADERS = libtetrabz.h libtetrabz.mod
lib_LTLIBRARIES = libtetrabz.la

libtetrabz_la_SOURCES = libtetrabz.F90

libtetrabz.o : libtetrabz.F90 libtetrabz_polcmplx_mod.o libtetrabz_fermigr_mod.o libtetrabz_polstat_mod.o libtetrabz_dbldelta_mod.o libtetrabz_dblstep_mod.o libtetrabz_dos_mod.o libtetrabz_occ_mod.o
libtetrabz_common.o : libtetrabz_common.F90
libtetrabz_dbldelta_mod.o : libtetrabz_dbldelta_mod.F90 libtetrabz_common.o
libtetrabz_dblstep_mod.o : libtetrabz_dblstep_mod.F90 libtetrabz_common.o
libtetrabz_dos_mod.o : libtetrabz_dos_mod.F90 libtetrabz_common.o
libtetrabz_fermigr_mod.o : libtetrabz_fermigr_mod.F90 libtetrabz_common.o
libtetrabz_occ_mod.o : libtetrabz_occ_mod.F90 libtetrabz_common.o
libtetrabz_polcmplx_mod.o : libtetrabz_polcmplx_mod.F90 libtetrabz_common.o
libtetrabz_polstat_mod.o : libtetrabz_polstat_mod.F90 libtetrabz_common.o
