Re-distribution of this program
===============================

Contain libtetrabz in your program
----------------------------------

libtetrabz is distributed with the :ref:`mitlicense`.
To summarize this, you can freely modify, copy and paste libtetrabz to any program
such as a private program (in the research group, co-workers, etc.),
open-source, free, and commercial software.
Also, you can freely choose the license to distribute your program.

Build libtetrabz without Autoconf
----------------------------------

In this package, libtetrabz is built with Autotools (Autoconf, Automake, Libtool).
If you do not want to use Autotools for your distributed program with libtetrabz's source,
you can use the following simple Makefile (please care about TAB).

.. code-block:: makefile

   F90 = gfortran
   FFLAGS = -fopenmp -O2 -g
   
   OBJS = \
   libtetrabz.o \
   libtetrabz_dbldelta_mod.o \
   libtetrabz_dblstep_mod.o \
   libtetrabz_dos_mod.o \
   libtetrabz_fermigr_mod.o \
   libtetrabz_occ_mod.o \
   libtetrabz_polcmplx_mod.o \
   libtetrabz_polstat_mod.o \
   libtetrabz_common.o

   .SUFFIXES :
   .SUFFIXES : .o .F90

   libtetrabz.a:$(OBJS)
        ar cr $@ $(OBJS)

   .F90.o:
         $(F90) $(FFLAGS) -c $<

   clean:
         rm -f *.a *.o *.mod

   libtetrabz.o:libtetrabz_polcmplx_mod.o
   libtetrabz.o:libtetrabz_fermigr_mod.o
   libtetrabz.o:libtetrabz_polstat_mod.o
   libtetrabz.o:libtetrabz_dbldelta_mod.o
   libtetrabz.o:libtetrabz_dblstep_mod.o
   libtetrabz.o:libtetrabz_dos_mod.o
   libtetrabz.o:libtetrabz_occ_mod.o
   libtetrabz_dbldelta_mod.o:libtetrabz_common.o
   libtetrabz_dblstep_mod.o:libtetrabz_common.o
   libtetrabz_dos_mod.o:libtetrabz_common.o
   libtetrabz_fermigr_mod.o:libtetrabz_common.o
   libtetrabz_occ_mod.o:libtetrabz_common.o
   libtetrabz_polcmplx_mod.o:libtetrabz_common.o
   libtetrabz_polstat_mod.o:libtetrabz_common.o

.. _mitlicense:

MIT License
-----------

| Copyright (c) 2014 Mitsuaki Kawamura
|
| Permission is hereby granted, free of charge, to any person obtaining a
| copy of this software and associated documentation files (the
| "Software"), to deal in the Software without restriction, including
| without limitation the rights to use, copy, modify, merge, publish,
| distribute, sublicense, and/or sell copies of the Software, and to
| permit persons to whom the Software is furnished to do so, subject to
| the following conditions:
| 
| The above copyright notice and this permission notice shall be included
| in all copies or substantial portions of the Software.
|
| THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
| OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
| MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
| IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
| CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
| TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
| SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
