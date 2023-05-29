プログラムの再配布
==================

自分のプログラムにlibtetrabzを含める
-------------------------------------

libtetrabzは下記の :ref:`mitlicense` に基づいて配布されている.
これはかいつまんで言うと,
個人的(研究室や共同研究者等のグループ)なプログラムであろうとも,
公開したり売ったりするプログラムであろうとも
自由にコピペしたり改変して良いし,
どのようなライセンスで配布しても構わない, と言うことである.

Autoconfを使わずにlibtetrabzをビルドする
-----------------------------------------

このパッケージではAutotools (Autoconf, Aitomake, Libtool)を使ってlibtetrabzをビルドしている.
もし再配布するソースコードにlibtetrabzを含めるときに,
Autoconfの使用に支障がある場合には, 以下の簡易版のMakefileを使うと良い (タブに注意).

.. code-block:: make

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
   libtetrabz_common.o \

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
   
MIT ライセンス
--------------

| Copyright (c) 2014 Mitsuaki Kawamura
| 
| 以下に定める条件に従い,
| 本ソフトウェアおよび関連文書のファイル（以下「ソフトウェア」）
| の複製を取得するすべての人に対し,
| ソフトウェアを無制限に扱うことを無償で許可します. これには,
| ソフトウェアの複製を使用, 複写, 変更, 結合, 掲載, 頒布, サブライセンス,
| および/または販売する権利,
| およびソフトウェアを提供する相手に同じことを許可する権利も無制限に含まれます.
| 
| 上記の著作権表示および本許諾表示を,
| ソフトウェアのすべての複製または重要な部分に記載するものとします.
| 
| ソフトウェアは「現状のまま」で, 明示であるか暗黙であるかを問わず,
| 何らの保証もなく提供されます. ここでいう保証とは, 商品性,
| 特定の目的への適合性, および権利非侵害についての保証も含みますが,
| それに限定されるものではありません. 作者または著作権者は, 契約行為,
| 不法行為, またはそれ以外であろうと, ソフトウェアに起因または関連し,
| あるいはソフトウェアの使用またはその他の扱いによって生じる一切の請求,
| 損害, その他の義務について何らの責任も負わないものとします.
