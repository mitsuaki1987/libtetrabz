/*
 Copyright (c) 2014 Mitsuaki Kawamura

 Permission is hereby granted, free of charge, to any person obtaining a
 copy of this software and associated documentation files (the
 "Software"), to deal in the Software without restriction, including
 without limitation the rights to use, copy, modify, merge, publish,
 distribute, sublicense, and/or sell copies of the Software, and to
 permit persons to whom the Software is furnished to do so, subject to
 the following conditions:
 
 The above copyright notice and this permission notice shall be included
 in all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void libtetrabz_initialize(
    int ng[3],
    double bvec[3][3],
    double wlsm[20][4],
    int **ikv
    )
    {
  /*
  define shortest diagonal line & define type of tetragonal
  */
  int ii, i0, i1, i2, i3, it, itype, ivvec0[4], divvec[4][4], ivvec[6][20][3], ikv0[3], nt;
  double bvec2[3][3], bvec3[4][3], bnorm[4];
  double wlsm1[4][4] = { {1440.0, 0.0, 30.0, 0.0},
                         {0.0, 1440.0, 0.0, 30.0},
                         {30.0, 0.0, 1440.0, 0.0},
                         {0.0, 30.0, 0.0, 1440.0} },
  wlsm2[4][4] = {{-38.0, 7.0, 17.0, -28.0},
                 {-28.0, -38.0, 7.0, 17.0},
                 {17.0, -28.0, -38.0, 7.0},
                 {7.0, 17.0, -28.0, -38.0}},
  wlsm3[4][4] = {{-56.0, 9.0, -46.0, 9.0},
                 {9.0, -56.0, 9.0, -46.0},
                 {-46.0, 9.0, -56.0, 9.0},
                 {9.0, -46.0, 9.0, -56.0}},
  wlsm4[4][4] = {{-38.0, -28.0, 17.0, 7.0},
                 {7.0, -38.0, -28.0, 17.0},
                 {17.0, 7.0, -38.0, -28.0},
                 {-28.0, 17.0, 7.0, -38.0}},
  wlsm5[4][4] = {{-18.0, -18.0, 12.0, -18.0},
                 {-18.0, -18.0, -18.0, 12.0},
                 {12.0, -18.0, -18.0, -18.0},
                 {-18.0, 12.0, -18.0, -18.0}};

  for (i1 = 0; i1 < 3; i1++)
    for (i2 = 0; i2 < 3; i2++)
      bvec2[i1][i2] = bvec[i1][i2] / (double) ng[i1];

  for (i1 = 0; i1 < 3; i1++) {
    bvec3[0][i1] = -bvec2[0][i1] + bvec2[1][i1] + bvec2[2][i1];
    bvec3[1][i1] = bvec2[0][i1] - bvec2[1][i1] + bvec2[2][i1];
    bvec3[2][i1] = bvec2[0][i1] + bvec2[1][i1] - bvec2[2][i1];
    bvec3[3][i1] = bvec2[0][i1] + bvec2[1][i1] + bvec2[2][i1];
  }
  /*
  length of delta bvec
  */
  for (i1 = 0; i1 < 4; i1++) {
    bnorm[i1] = 0.0;
    for (i2 = 0; i2 < 3; i2++)
      bnorm[i1] += bvec3[i1][i2] * bvec3[i1][i2];
  }
  itype = 0;
  for (i1 = 1; i1 < 4; i1++)
    if (bnorm[i1] < bnorm[itype]) itype = i1;
  /*
  start & last
  */
  for (i0 = 0; i0 < 4; i0++) {
    ivvec0[i0] = 0;
    for (i1 = 0; i1 < 4; i1++)divvec[i0][i1] = 0;
    divvec[i0][i0] = 1;
  }
  ivvec0[itype] = 1;
  divvec[itype][itype] = -1;
  /*
  Corners of tetrahedron
  */
  it = 0;
  for (i0 = 0; i0 < 3; i0++) {
    for (i1 = 0; i1 < 3; i1++) {
      if (i1 == i0) continue;
      for (i2 = 0; i2 < 3; i2++) {
        if (i2 == i1 || i2 == i0) continue;

        for (i3 = 0; i3 < 3; i3++) {
          ivvec[it][0][i3] = ivvec0[i3];
          ivvec[it][1][i3] = ivvec[it][0][i3] + divvec[i0][i3];
          ivvec[it][2][i3] = ivvec[it][1][i3] + divvec[i1][i3];
          ivvec[it][3][i3] = ivvec[it][2][i3] + divvec[i2][i3];
        }

        it += 1;
      }
    }
  }
  /*
  Additional points
  */
  for (i1 = 0; i1 < 6; i1++) {
    for (i2 = 0; i2 < 3; i2++) {
      ivvec[i1][4][i2] = 2 * ivvec[i1][0][i2] - ivvec[i1][1][i2];
      ivvec[i1][5][i2] = 2 * ivvec[i1][1][i2] - ivvec[i1][2][i2];
      ivvec[i1][6][i2] = 2 * ivvec[i1][2][i2] - ivvec[i1][3][i2];
      ivvec[i1][7][i2] = 2 * ivvec[i1][3][i2] - ivvec[i1][0][i2];

      ivvec[i1][8][i2] = 2 * ivvec[i1][0][i2] - ivvec[i1][2][i2];
      ivvec[i1][9][i2] = 2 * ivvec[i1][1][i2] - ivvec[i1][3][i2];
      ivvec[i1][10][i2] = 2 * ivvec[i1][2][i2] - ivvec[i1][0][i2];
      ivvec[i1][11][i2] = 2 * ivvec[i1][3][i2] - ivvec[i1][1][i2];

      ivvec[i1][12][i2] = 2 * ivvec[i1][0][i2] - ivvec[i1][3][i2];
      ivvec[i1][13][i2] = 2 * ivvec[i1][1][i2] - ivvec[i1][0][i2];
      ivvec[i1][14][i2] = 2 * ivvec[i1][2][i2] - ivvec[i1][1][i2];
      ivvec[i1][15][i2] = 2 * ivvec[i1][3][i2] - ivvec[i1][2][i2];

      ivvec[i1][16][i2] = ivvec[i1][3][i2] - ivvec[i1][0][i2] + ivvec[i1][1][i2];
      ivvec[i1][17][i2] = ivvec[i1][0][i2] - ivvec[i1][1][i2] + ivvec[i1][2][i2];
      ivvec[i1][18][i2] = ivvec[i1][1][i2] - ivvec[i1][2][i2] + ivvec[i1][3][i2];
      ivvec[i1][19][i2] = ivvec[i1][2][i2] - ivvec[i1][3][i2] + ivvec[i1][0][i2];
    }
  }

  for (i1 = 0; i1 < 4; i1++) {
    for (i2 = 0; i2 < 4; i2++) {
      wlsm[i2][i1] = wlsm1[i1][i2] /= 1260.0;
      wlsm[i2 + 4][i1] = wlsm2[i1][i2] /= 1260.0;
      wlsm[i2 + 8][i1] = wlsm3[i1][i2] /= 1260.0;
      wlsm[i2 + 12][i1] = wlsm4[i1][i2] /= 1260.0;
      wlsm[i2 + 16][i1] = wlsm5[i1][i2] /= 1260.0;
    }
  }
  /*
  k-index for energy
  */
  nt = 0;
  for (i2 = 0; i2 < ng[2]; i2++) {
    for (i1 = 0; i1 < ng[1]; i1++) {
      for (i0 = 0; i0 < ng[0]; i0++) {

        for (it = 0; it < 6; it++) {

          for (ii = 0; ii < 20; ii++) {
            ikv0[0] = (i0 + ivvec[it][ii][0]) % ng[0];
            ikv0[1] = (i1 + ivvec[it][ii][1]) % ng[1];
            ikv0[2] = (i2 + ivvec[it][ii][2]) % ng[2];
            for (i3 = 0; i3 < 3; i3++) if (ikv0[i3] < 0) ikv0[i3] += ng[i3];
            ikv[nt][ii] = ikv0[2] + ng[2] * ikv0[1] + ng[2] * ng[1] * ikv0[0];
          }
          nt += 1;
        }
      }
    }
  }
}
/*
Cut small tetrahedron A1
*/
void libtetrabz_tsmall_a1(
    double *e,
    double e0,
    double *v,
    double tsmall[4][4]
) {
  double a10, a20, a30;
  a10 = (e0 - e[0]) / (e[1] - e[0]);
  a20 = (e0 - e[0]) / (e[2] - e[0]);
  a30 = (e0 - e[0]) / (e[3] - e[0]);

  *v = a10 * a20 * a30;

  tsmall[0][0] = 1.0;
  tsmall[0][1] = 1.0 - a10;
  tsmall[0][2] = 1.0 - a20;
  tsmall[0][3] = 1.0 - a30;
  tsmall[1][0] = 0.0;
  tsmall[1][1] = a10;
  tsmall[1][2] = 0.0;
  tsmall[1][3] = 0.0;
  tsmall[2][0] = 0.0;
  tsmall[2][1] = 0.0;
  tsmall[2][2] = a20;
  tsmall[2][3] = 0.0;
  tsmall[3][0] = 0.0;
  tsmall[3][1] = 0.0;
  tsmall[3][2] = 0.0;
  tsmall[3][3] = a30;
}
/*
Cut small tetrahedron B1
*/
void libtetrabz_tsmall_b1(
  double *e,
  double e0,
  double *v,
  double tsmall[4][4]
)
{
  double a13, a20, a30;
  a13 = (e0 - e[3]) / (e[1] - e[3]);
  a20 = (e0 - e[0]) / (e[2] - e[0]);
  a30 = (e0 - e[0]) / (e[3] - e[0]);

  *v = a20 * a30 * a13;

  tsmall[0][0] = 1.0;
  tsmall[0][1] = 1.0 - a20;
  tsmall[0][2] = 1.0 - a30;
  tsmall[0][3] = 0.0;
  tsmall[1][0] = 0.0;
  tsmall[1][1] = 0.0;
  tsmall[1][2] = 0.0;
  tsmall[1][3] = a13;
  tsmall[2][0] = 0.0;
  tsmall[2][1] = a20;
  tsmall[2][2] = 0.0;
  tsmall[2][3] = 0.0;
  tsmall[3][0] = 0.0;
  tsmall[3][1] = 0.0;
  tsmall[3][2] = a30;
  tsmall[3][3] = 1.0 - a13;
}
/*
Cut small tetrahedron B2
*/
void libtetrabz_tsmall_b2(
  double *e,
  double e0,
  double *v,
  double tsmall[4][4]
)
{
  double a21, a31;
  a21 = (e0 - e[1]) / (e[2] - e[1]);
  a31 = (e0 - e[1]) / (e[3] - e[1]);

  *v = a21 * a31;

  tsmall[0][0] = 1.0;
  tsmall[0][1] = 0.0;
  tsmall[0][2] = 0.0;
  tsmall[0][3] = 0.0;
  tsmall[1][0] = 0.0;
  tsmall[1][1] = 1.0;
  tsmall[1][2] = 1.0 - a21;
  tsmall[1][3] = 1.0 - a31;
  tsmall[2][0] = 0.0;
  tsmall[2][1] = 0.0;
  tsmall[2][2] = a21;
  tsmall[2][3] = 0.0;
  tsmall[3][0] = 0.0;
  tsmall[3][1] = 0.0;
  tsmall[3][2] = 0.0;
  tsmall[3][3] = a31;
}
/*
Cut small tetrahedron B3
*/
void libtetrabz_tsmall_b3(
  double *e,
  double e0,
  double *v,
  double tsmall[4][4]
)
{
  double a12, a20, a31;
  a12 = (e0 - e[2]) / (e[1] - e[2]);
  a20 = (e0 - e[0]) / (e[2] - e[0]);
  a31 = (e0 - e[1]) / (e[3] - e[1]);

  *v = a12 * a20 * a31;

  tsmall[0][0] = 1.0;
  tsmall[0][1] = 1.0 - a20;
  tsmall[0][2] = 0.0;
  tsmall[0][3] = 0.0;
  tsmall[1][0] = 0.0;
  tsmall[1][1] = 0.0;
  tsmall[1][2] = a12;
  tsmall[1][3] = 1.0 - a31;
  tsmall[2][0] = 0.0;
  tsmall[2][1] = a20;
  tsmall[2][2] = 1.0 - a12;
  tsmall[2][3] = 0.0;
  tsmall[3][0] = 0.0;
  tsmall[3][1] = 0.0;
  tsmall[3][2] = 0.0;
  tsmall[3][3] = a31;
}
/*
Cut small tetrahedron C1
*/
void libtetrabz_tsmall_c1(
  double *e,
  double e0,
  double *v,
  double tsmall[4][4]
)
{
  double a32;
  a32 = (e0 - e[2]) / (e[3] - e[2]);

  *v = a32;

  tsmall[0][0] = 1.0;
  tsmall[0][1] = 0.0;
  tsmall[0][2] = 0.0;
  tsmall[0][3] = 0.0;
  tsmall[1][0] = 0.0;
  tsmall[1][1] = 1.0;
  tsmall[1][2] = 0.0;
  tsmall[1][3] = 0.0;
  tsmall[2][0] = 0.0;
  tsmall[2][1] = 0.0;
  tsmall[2][2] = 1.0;
  tsmall[2][3] = 1.0 - a32;
  tsmall[3][0] = 0.0;
  tsmall[3][1] = 0.0;
  tsmall[3][2] = 0.0;
  tsmall[3][3] = a32;
}
/*
Cut small tetrahedron C2
*/
void libtetrabz_tsmall_c2(
  double *e,
  double e0,
  double *v,
  double tsmall[4][4]
)
{
  double a23, a31;
  a23 = (e0 - e[3]) / (e[2] - e[3]);
  a31 = (e0 - e[1]) / (e[3] - e[1]);

  *v = a23 * a31;

  tsmall[0][0] = 1.0;
  tsmall[0][1] = 0.0;
  tsmall[0][2] = 0.0;
  tsmall[0][3] = 0.0;
  tsmall[1][0] = 0.0;
  tsmall[1][1] = 1.0;
  tsmall[1][2] = 1.0 - a31;
  tsmall[1][3] = 0.0;
  tsmall[2][0] = 0.0;
  tsmall[2][1] = 0.0;
  tsmall[2][2] = 0.0;
  tsmall[2][3] = a23;
  tsmall[3][0] = 0.0;
  tsmall[3][1] = 0.0;
  tsmall[3][2] = a31;
  tsmall[3][3] = 1.0 - a23;
}
/*
Cut small tetrahedron C3
*/
void libtetrabz_tsmall_c3(
  double *e,
  double e0,
  double *v,
  double tsmall[4][4]
)
{
  double a23, a13, a30;
  a23 = (e0 - e[3]) / (e[2] - e[3]);
  a13 = (e0 - e[3]) / (e[1] - e[3]);
  a30 = (e0 - e[0]) / (e[3] - e[0]);

  *v = a23 * a13 * a30;

  tsmall[0][0] = 1.0;
  tsmall[0][1] = 1.0 - a30;
  tsmall[0][2] = 0.0;
  tsmall[0][3] = 0.0;
  tsmall[1][0] = 0.0;
  tsmall[1][1] = 0.0;
  tsmall[1][2] = a13;
  tsmall[1][3] = 0.0;
  tsmall[2][0] = 0.0;
  tsmall[2][1] = 0.0;
  tsmall[2][2] = 0.0;
  tsmall[2][3] = a23;
  tsmall[3][0] = 0.0;
  tsmall[3][1] = a30;
  tsmall[3][2] = 1.0 - a13;
  tsmall[3][3] = 1.0 - a23;
}
/*
Cut triangle A1
*/
void libtetrabz_triangle_a1(
  double *e,
  double e0,
  double *v,
  double tsmall[4][3]
)
{
  double a10, a20, a30;
  a10 = (e0 - e[0]) / (e[1] - e[0]);
  a20 = (e0 - e[0]) / (e[2] - e[0]);
  a30 = (e0 - e[0]) / (e[3] - e[0]);

  *v = 3.0 * a10 * a20 / (e[3] - e[0]);

  tsmall[0][0] = 1.0 - a10;
  tsmall[0][1] = 1.0 - a20;
  tsmall[0][2] = 1.0 - a30;
  tsmall[1][0] = a10;
  tsmall[1][1] = 0.0;
  tsmall[1][2] = 0.0;
  tsmall[2][0] = 0.0;
  tsmall[2][1] = a20;
  tsmall[2][2] = 0.0;
  tsmall[3][0] = 0.0;
  tsmall[3][1] = 0.0;
  tsmall[3][2] = a30;
}
/*
Cut triangle B1
*/
void libtetrabz_triangle_b1(
  double *e,
  double e0,
  double *v,
  double tsmall[4][3]
) {
  double a30, a13, a20;
  a30 = (e0 - e[0]) / (e[3] - e[0]);
  a13 = (e0 - e[3]) / (e[1] - e[3]);
  a20 = (e0 - e[0]) / (e[2] - e[0]);

  *v = 3.0 * a30 * a13 / (e[2] - e[0]);

  tsmall[0][0] = 1.0 - a20;
  tsmall[0][1] = 1.0 - a30;
  tsmall[0][2] = 0.0;
  tsmall[1][0] = 0.0;
  tsmall[1][1] = 0.0;
  tsmall[1][2] = a13;
  tsmall[2][0] = a20;
  tsmall[2][1] = 0.0;
  tsmall[2][2] = 0.0;
  tsmall[3][0] = 0.0;
  tsmall[3][1] = a30;
  tsmall[3][2] = 1.0 - a13;
}
/*
Cut triangle B2
*/
void libtetrabz_triangle_b2(
    double *e,
    double e0,
    double *v,
    double tsmall[4][3]
)
{
  double a12, a31, a20;
  a12 = (e0 - e[2]) / (e[1] - e[2]);
  a31 = (e0 - e[1]) / (e[3] - e[1]);
  a20 = (e0 - e[0]) / (e[2] - e[0]);

  *v = 3.0 * a12 * a31 / (e[2] - e[0]);

  tsmall[0][0] = 1.0 - a20;
  tsmall[0][1] = 0.0;
  tsmall[0][2] = 0.0;
  tsmall[1][0] = 0.0;
  tsmall[1][1] = a12;
  tsmall[1][2] = 1.0 - a31;
  tsmall[2][0] = a20;
  tsmall[2][1] = 1.0 - a12;
  tsmall[2][2] = 0.0;
  tsmall[3][0] = 0.0;
  tsmall[3][1] = 0.0;
  tsmall[3][2]= a31;
}
/*
Cut triangle C1
*/
void libtetrabz_triangle_c1(
  double *e,
  double e0,
  double *v,
  double tsmall[4][3]
)
{
  double a03, a13, a23;
  a03 = (e0 - e[3]) / (e[0] - e[3]);
  a13 = (e0 - e[3]) / (e[1] - e[3]);
  a23 = (e0 - e[3]) / (e[2] - e[3]);

  *v = 3.0 * a03 * a13 / (e[3] - e[2]);

  tsmall[0][0] = a03;
  tsmall[0][1] = 0.0;
  tsmall[0][2] = 0.0;
  tsmall[1][0] = 0.0;
  tsmall[1][1] = a13;
  tsmall[1][2] = 0.0;
  tsmall[2][0] = 0.0;
  tsmall[2][1] = 0.0;
  tsmall[2][2] = a23;
  tsmall[3][0] = 1.0 - a03;
  tsmall[3][1] = 1.0 - a13;
  tsmall[3][2] = 1.0 - a23;
}
/*
Sort eigenvalues
*/
void eig_sort(
    int n, //!< [in] the number of components
    double *key, //!< [in] Variables to be sorted [n].
    int *swap //!< [out] Order of index (sorted)
)
{
  int i, j, k, min_loc;
  double min_val;

  for (i = 0; i < n; ++i) swap[i] = i;

  for (i = 0; i < n - 1; ++i) {
    min_val = key[i];
    min_loc = i;
    for (j = i + 1; j < n; ++j) {
      if (min_val > key[j]) {
        min_val = key[j];
        min_loc = j;
      }
    }
    if (key[i] > min_val) {
      /*
       Swap
      */
      key[min_loc] = key[i];
      key[i] = min_val;

      k = swap[min_loc];
      swap[min_loc] = swap[i];
      swap[i] = k;
    }
  }/*for (i = 0; i < n - 1; ++i)*/
}/*eig_sort*/
/*
2nd step of tetrahedron method.
*/
static void libtetrabz_dbldelta2(
    int nb,
    double **ej,
    double **w
) {
  int i3, ib, indx[3];
  double a10, a20, a02, a12, v, e[3], e_abs;

  for (ib = 0; ib < nb; ib++) {

    e_abs = 0.0;
    for (i3 = 0; i3 < 3; i3++) {
      e[i3] = ej[ib][i3];
      if (e_abs < fabs(e[i3])) e_abs = fabs(e[i3]);
    }
    eig_sort(3, e, indx);

    if (e_abs < 1.0e-10) {
      printf("Nesting ##\n");
    }

    if ((e[0] < 0.0 && 0.0 <= e[1]) || (e[0] <= 0.0 && 0.0 < e[1])) {

      a10 = (0.0 - e[0]) / (e[1] - e[0]);
      a20 = (0.0 - e[0]) / (e[2] - e[0]);

      v = a10 / (e[2] - e[0]);

      w[ib][indx[0]] = v * (2.0 - a10 - a20);
      w[ib][indx[1]] = v * a10;
      w[ib][indx[2]] = v * a20;
    }
    else if ((e[1] <= 0.0 && 0.0 < e[2]) || (e[1] < 0.0 && 0.0 <= e[2])) {

      a02 = (0.0 - e[2]) / (e[0] - e[2]);
      a12 = (0.0 - e[2]) / (e[1] - e[2]);

      v = a12 / (e[2] - e[0]);

      w[ib][indx[0]] = v * a02;
      w[ib][indx[1]] = v * a12;
      w[ib][indx[2]] = v * (2.0 - a02 - a12);
    }
    else {
      for (i3 = 0; i3 < 3; i3++)
        w[ib][i3] = 0.0;
    }
  }
}
/*
Main SUBROUTINE for Delta(E1) * Delta(E2)
*/
static PyObject* dbldelta_c(PyObject* self, PyObject* args)
{
  int it, ik, ib, i20, i4, j3, jb, ** ikv, indx[4], ierr, ng[3], nk, nb;
  double wlsm[20][4], ** ei1, ** ej1, ** ej2, e[4], *** w1, **w2, v, tsmall[4][3], thr, 
    bvec[3][3], ** eig1, ** eig2, *** wght;
  PyObject* eig1_po, * eig2_po, * wght_po;
  /*
   Read input from python object
   */
  if (!PyArg_ParseTuple(args, "iiiiidddddddddOO",
    &ng[0], &ng[1], &ng[2], &nk, &nb,
    &bvec[0][0], &bvec[0][1], &bvec[0][2], &bvec[1][0], &bvec[1][1], &bvec[1][2], &bvec[2][0], &bvec[2][1], &bvec[2][2],
    &eig1_po, &eig2_po))
    return NULL;
  /*
  convert python list object to array
  */
  eig1 = (double**)malloc(nk * sizeof(double*));
  eig1[0] = (double*)malloc(nk * nb * sizeof(double));
  eig2 = (double**)malloc(nk * sizeof(double*));
  eig2[0] = (double*)malloc(nk * nb * sizeof(double));
  wght = (double***)malloc(nk * sizeof(double**));
  wght[0] = (double**)malloc(nk * nb * sizeof(double*));
  wght[0][0] = (double*)malloc(nk * nb * nb * sizeof(double));
  for (ik = 0; ik < nk; ik++) {
    eig1[ik] = eig1[0] + ik * nb;
    eig2[ik] = eig2[0] + ik * nb;
    wght[ik] = wght[0] + ik * nb;
    for (ib = 0; ib < nb; ib++) {
      wght[ik][ib] = wght[0][0] + ik * nb * nb + ib * nb;
    }
  }

  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++) {
      eig1[ik][ib] = PyFloat_AsDouble(PyList_GetItem(eig1_po, ik * nb + ib));
      eig2[ik][ib] = PyFloat_AsDouble(PyList_GetItem(eig2_po, ik * nb + ib));
    }
  /*
  Start main calculation
  */
  thr = 1.0e-10;

  ikv = (int**)malloc(6 * nk * sizeof(int*));
  ikv[0] = (int*)malloc(6 * nk * 20 * sizeof(int));
  for (ik = 0; ik < 6 * nk; ik++) {
    ikv[ik] = ikv[0] + ik * 20;
  }

  ei1 = (double**)malloc(4 * sizeof(double*));
  ej1 = (double**)malloc(4 * sizeof(double*));
  ei1[0] = (double*)malloc(4 * nb * sizeof(double));
  ej1[0] = (double*)malloc(4 * nb * sizeof(double));
  for (i4 = 0; i4 < 4; i4++) {
    ei1[i4] = ei1[0] + i4 * nb;
    ej1[i4] = ej1[0] + i4 * nb;
  }

  w1 = (double***)malloc(nb * sizeof(double**));
  w1[0] = (double**)malloc(nb * 4 * sizeof(double*));
  w1[0][0] = (double*)malloc(nb * 4 * nb * sizeof(double));
  for (ib = 0; ib < nb; ib++) {
    w1[ib] = w1[0] + ib * 4;
    for (i4 = 0; i4 < 4; i4++) {
      w1[ib][i4] = w1[0][0] + ib * 4 * nb + i4 * nb;
    }
  }

  ej2 = (double**)malloc(nb * sizeof(double*));
  ej2[0] = (double*)malloc(nb * 3 * sizeof(double));
  for (ib = 0; ib < nb; ib++) {
    ej2[ib] = ej2[0] + ib * 3;
  }

  w2 = (double**)malloc(nb * sizeof(double*));
  w2[0] = (double*)malloc(nb * 3 * sizeof(double));
  for (ib = 0; ib < nb; ib++) {
    w2[ib] = w2[0] + ib * 3;
  }

  libtetrabz_initialize(ng, bvec, wlsm, ikv);

  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      for (jb = 0; jb < nb; jb++)
        wght[ik][ib][jb] = 0.0;

  for (it = 0; it < 6 * nk; it++) {

    for (i4 = 0; i4 < 4; i4++)
      for (ib = 0; ib < nb; ib++) {
        ei1[i4][ib] = 0.0;
        ej1[i4][ib] = 0.0;
      }
    for (i20 = 0; i20 < 20; i20++) {
      for (i4 = 0; i4 < 4; i4++) {
        for (ib = 0; ib < nb; ib++) {
          ei1[i4][ib] += eig1[ikv[it][i20]][ib] * wlsm[i20][i4];
          ej1[i4][ib] += eig2[ikv[it][i20]][ib] * wlsm[i20][i4];
        }
      }
    }

    for (ib = 0; ib < nb; ib++) {

      for (i4 = 0; i4 < 4; i4++)
        for (jb = 0; jb < nb; jb++)
          w1[ib][i4][jb] = 0.0;

      for (i4 = 0; i4 < 4; i4++) e[i4] = ei1[i4][ib];
      eig_sort(4, e, indx);

      if (e[0] < 0.0 && 0.0 <= e[1]) {

        libtetrabz_triangle_a1(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (jb = 0; jb < nb; jb++)
            for (j3 = 0; j3 < 3; j3++) ej2[jb][j3] = 0.0;
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j3 = 0; j3 < 3; j3++)
                ej2[jb][j3] += ej1[indx[i4]][jb] * tsmall[i4][j3];
          libtetrabz_dbldelta2(nb, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j3 = 0; j3 < 3; j3++)
                w1[ib][indx[i4]][jb] += v * tsmall[i4][j3] * w2[jb][j3];
        }
      }
      else if (e[1] < 0.0 && 0.0 <= e[2]) {

        libtetrabz_triangle_b1(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (jb = 0; jb < nb; jb++)
            for (j3 = 0; j3 < 3; j3++) ej2[jb][j3] = 0.0;
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j3 = 0; j3 < 3; j3++)
                ej2[jb][j3] += ej1[indx[i4]][jb] * tsmall[i4][j3];
          libtetrabz_dbldelta2(nb, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j3 = 0; j3 < 3; j3++)
                w1[ib][indx[i4]][jb] += v * tsmall[i4][j3] * w2[jb][j3];
        }

        libtetrabz_triangle_b2(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (jb = 0; jb < nb; jb++)
            for (j3 = 0; j3 < 3; j3++) ej2[jb][j3] = 0.0;
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j3 = 0; j3 < 3; j3++)
                ej2[jb][j3] += ej1[indx[i4]][jb] * tsmall[i4][j3];
          libtetrabz_dbldelta2(nb, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j3 = 0; j3 < 3; j3++)
                w1[ib][indx[i4]][jb] += v * tsmall[i4][j3] * w2[jb][j3];
        }
      }
      else if (e[2] < 0.0 && 0.0 < e[3]) {

        libtetrabz_triangle_c1(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (jb = 0; jb < nb; jb++)
            for (j3 = 0; j3 < 3; j3++) ej2[jb][j3] = 0.0;
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j3 = 0; j3 < 3; j3++)
                ej2[jb][j3] += ej1[indx[i4]][jb] * tsmall[i4][j3];
          libtetrabz_dbldelta2(nb, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j3 = 0; j3 < 3; j3++)
                w1[ib][indx[i4]][jb] += v * tsmall[i4][j3] * w2[jb][j3];
        }
      }
      else {
        continue;
      }
    }
    for (i20 = 0; i20 < 20; i20++)
      for (ib = 0; ib < nb; ib++)
        for (i4 = 0; i4 < 4; i4++)
          for (jb = 0; jb < nb; jb++)
            wght[ikv[it][i20]][ib][jb] += wlsm[i20][i4] * w1[ib][i4][jb];
  }
  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      for (jb = 0; jb < nb; jb++)
        wght[ik][ib][jb] /= (6.0 * (double)nk);

  free(ikv[0]);
  free(ikv);
  free(ei1[0]);
  free(ei1);
  free(ej1[0]);
  free(ej1);
  free(ej2[0]);
  free(ej2);
  free(w1[0][0]);
  free(w1[0]);
  free(w1); 
  free(w2[0]);
  free(w2);
  /*
  Convert weight to python list object
  */
  wght_po = PyList_New(nk * nb * nb);
  ierr = 0;
  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      for (jb = 0; jb < nb; jb++)
        ierr = PyList_SetItem(wght_po, ik * nb * nb + ib * nb + jb, PyFloat_FromDouble(wght[ik][ib][jb]));
  if (ierr != 0)printf("Error in PyList_SetItem\n");

  free(eig1[0]);
  free(eig1);
  free(eig2[0]);
  free(eig2);
  free(wght[0][0]);
  free(wght[0]);
  free(wght);

  return wght_po;
}
/*
Tetrahedron method for theta( - de)
*/
static void libtetrabz_dblstep2(
    int nb,
    double ei[4],
    double **ej,
    double **w1
) {
  int i4, j4, ib, indx[4];
  double v, thr, e[4], tsmall[4][4];

  thr = 1.0e-8;

  for (ib = 0; ib < nb; ib++) {

    for (i4 = 0; i4 < 4; i4++)
      w1[ib][i4] = 0.0;

    for (i4 = 0; i4 < 4; i4++) e[i4] = - ei[i4] + ej[ib][i4];
    eig_sort(4, e, indx);

    if (fabs(e[0]) < thr && fabs(e[3]) < thr) {
      /*
      Theta(0) = 0.5
      */
      v = 0.5;
      for (i4 = 0; i4 < 4; i4++)
        w1[ib][i4] += v * 0.25;
    }
    else if ((e[0] <= 0.0 && 0.0 < e[1]) || (e[0] < 0.0 && 0.0 <= e[1])) {
      libtetrabz_tsmall_a1(e, 0.0, &v, tsmall);
      for (i4 = 0; i4 < 4; i4++)
        for (j4 = 0; j4 < 4; j4++)
          w1[ib][indx[i4]] += v * tsmall[i4][j4] * 0.25;
    }
    else if ((e[1] <= 0.0 && 0.0 < e[2]) || (e[1] < 0.0 && 0.0 <= e[2])) {

      libtetrabz_tsmall_b1(e, 0.0, &v, tsmall);
      for (i4 = 0; i4 < 4; i4++)
        for (j4 = 0; j4 < 4; j4++)
          w1[ib][indx[i4]] += v * tsmall[i4][j4] * 0.25;

      libtetrabz_tsmall_b2(e, 0.0, &v, tsmall);
      for (i4 = 0; i4 < 4; i4++)
        for (j4 = 0; j4 < 4; j4++)
          w1[ib][indx[i4]] += v * tsmall[i4][j4] * 0.25;

      libtetrabz_tsmall_b3(e, 0.0, &v, tsmall);
      for (i4 = 0; i4 < 4; i4++)
        for (j4 = 0; j4 < 4; j4++)
          w1[ib][indx[i4]] += v * tsmall[i4][j4] * 0.25;
    }
    else if ((e[2] <= 0.0 && 0.0 < e[3]) || (e[2] < 0.0 && 0.0 <= e[3])) {

      libtetrabz_tsmall_c1(e, 0.0, &v, tsmall);
      for (i4 = 0; i4 < 4; i4++)
        for (j4 = 0; j4 < 4; j4++)
          w1[ib][indx[i4]] += v * tsmall[i4][j4] * 0.25;

      libtetrabz_tsmall_c2(e, 0.0, &v, tsmall);
      for (i4 = 0; i4 < 4; i4++)
        for (j4 = 0; j4 < 4; j4++)
          w1[ib][indx[i4]] += v * tsmall[i4][j4] * 0.25;

      libtetrabz_tsmall_c3(e, 0.0, &v, tsmall);
      for (i4 = 0; i4 < 4; i4++)
        for (j4 = 0; j4 < 4; j4++)
          w1[ib][indx[i4]] += v * tsmall[i4][j4] * 0.25;
    }
    else if (e[3] <= 0.0) {
      for (i4 = 0; i4 < 4; i4++)
        w1[ib][i4] += 0.25;
    }
  }
}
/*
Main SUBROUTINE for Theta(- E1) * Theta(E1 - E2)
*/
static PyObject* dblstep_c(PyObject* self, PyObject* args)
{
  int it, ik, ib, jb, i20, i4, j4, **ikv, indx[4], ierr, ng[3], nk, nb;
  double wlsm[20][4], **ei1, **ej1, ei2[4], ** ej2, e[4], *** w1, ** w2, v, tsmall[4][4], thr,
    bvec[3][3], ** eig1, ** eig2, *** wght;
  PyObject* eig1_po, * eig2_po, * wght_po;

  /*
   Read input from python object
   */
  if (!PyArg_ParseTuple(args, "iiiiidddddddddOO",
    &ng[0], &ng[1], &ng[2], &nk, &nb,
    &bvec[0][0], &bvec[0][1], &bvec[0][2], &bvec[1][0], &bvec[1][1], &bvec[1][2], &bvec[2][0], &bvec[2][1], &bvec[2][2],
    &eig1_po, &eig2_po))
    return NULL;
  /*
  convert python list object to array
  */
  eig1 = (double**)malloc(nk * sizeof(double*));
  eig1[0] = (double*)malloc(nk * nb * sizeof(double));
  eig2 = (double**)malloc(nk * sizeof(double*));
  eig2[0] = (double*)malloc(nk * nb * sizeof(double));
  wght = (double***)malloc(nk * sizeof(double**));
  wght[0] = (double**)malloc(nk * nb * sizeof(double*));
  wght[0][0] = (double*)malloc(nk * nb * nb * sizeof(double));
  for (ik = 0; ik < nk; ik++) {
    eig1[ik] = eig1[0] + ik * nb;
    eig2[ik] = eig2[0] + ik * nb;
    wght[ik] = wght[0] + ik * nb;
    for (ib = 0; ib < nb; ib++) {
      wght[ik][ib] = wght[0][0] + ik * nb * nb + ib * nb;
    }
  }

  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++) {
      eig1[ik][ib] = PyFloat_AsDouble(PyList_GetItem(eig1_po, ik * nb + ib));
      eig2[ik][ib] = PyFloat_AsDouble(PyList_GetItem(eig2_po, ik * nb + ib));
    }
  /*
  Start main calculation
  */
  ikv = (int**)malloc(6 * nk * sizeof(int*));
  ikv[0] = (int*)malloc(6 * nk * 20 * sizeof(int));
  for (ik = 0; ik < 6 * nk; ik++) {
    ikv[ik] = ikv[0] + ik * 20;
  }

  ei1 = (double**)malloc(4 * sizeof(double*));
  ej1 = (double**)malloc(4 * sizeof(double*));
  ei1[0] = (double*)malloc(4 * nb * sizeof(double));
  ej1[0] = (double*)malloc(4 * nb * sizeof(double));
  for (i4 = 0; i4 < 4; i4++) {
    ei1[i4] = ei1[0] + i4 * nb;
    ej1[i4] = ej1[0] + i4 * nb;
  }

  ej2 = (double**)malloc(nb * sizeof(double*));
  ej2[0] = (double*)malloc(nb * 4 * sizeof(double));
  for (ib = 0; ib < nb; ib++) {
    ej2[ib] = ej2[0] + ib * 4;
  }

  w1 = (double***)malloc(nb * sizeof(double**));
  w1[0] = (double**)malloc(nb * 4 * sizeof(double*));
  w1[0][0] = (double*)malloc(nb * 4 * nb * sizeof(double));
  for (ib = 0; ib < nb; ib++) {
    w1[ib] = w1[0] + ib * 4;
    for (i4 = 0; i4 < 4; i4++) {
      w1[ib][i4] = w1[0][0] + ib * 4 * nb + i4 * nb;
    }
  }

  w2 = (double**)malloc(nb * sizeof(double*));
  w2[0] = (double*)malloc(nb * 4 * sizeof(double));
  for (ib = 0; ib < nb; ib++) {
    w2[ib] = w2[0] + ib * 4;
  }

  libtetrabz_initialize(ng, bvec, wlsm, ikv);

  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      for (jb = 0; jb < nb; jb++)
        wght[ik][ib][jb] = 0.0;

  thr = 1.0e-8;

  for(it = 0; it < 6*nk; it++) {

    for (i4 = 0; i4 < 4; i4++)
      for (ib = 0; ib < nb; ib++) {
        ei1[i4][ib] = 0.0;
        ej1[i4][ib] = 0.0;
      }
    for (i20 = 0; i20 < 20; i20++) {
      for (i4 = 0; i4 < 4; i4++) {
        for (ib = 0; ib < nb; ib++) {
          ei1[i4][ib] += eig1[ikv[it][i20]][ib] * wlsm[i20][i4];
          ej1[i4][ib] += eig2[ikv[it][i20]][ib] * wlsm[i20][i4];
        }
      }
    }
    
    for (ib = 0; ib < nb; ib++) {

      for (i4 = 0; i4 < 4; i4++)
        for (jb = 0; jb < nb; jb++)
          w1[ib][i4][jb] = 0.0;

      for (i4 = 0; i4 < 4; i4++) e[i4] = ei1[i4][ib];
      eig_sort(4, e, indx);

      if (e[0] <= 0.0 && 0.0 < e[1]) {

        libtetrabz_tsmall_a1(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (j4 = 0; j4 < 4; j4++) {
            ei2[j4] = 0.0;
            for (jb = 0; jb < nb; jb++) ej2[jb][j4] = 0.0;
          }
          for (i4 = 0; i4 < 4; i4++)
            for (j4 = 0; j4 < 4; j4++) {
              ei2[j4] += ei1[indx[i4]][ib] * tsmall[i4][j4];
              for (jb = 0; jb < nb; jb++)
                ej2[jb][j4] += ej1[indx[i4]][jb] * tsmall[i4][j4];
            }
          libtetrabz_dblstep2(nb, ei2, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j4 = 0; j4 < 4; j4++)
                w1[ib][indx[i4]][jb] += v * tsmall[i4][j4] * w2[jb][j4];
        }
      }
      else if (e[1] <= 0.0 && 0.0 < e[2]) {

        libtetrabz_tsmall_b1(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (j4 = 0; j4 < 4; j4++) {
            ei2[j4] = 0.0;
            for (jb = 0; jb < nb; jb++) ej2[jb][j4] = 0.0;
          }
          for (i4 = 0; i4 < 4; i4++)
            for (j4 = 0; j4 < 4; j4++) {
              ei2[j4] += ei1[indx[i4]][ib] * tsmall[i4][j4];
              for (jb = 0; jb < nb; jb++)
                ej2[jb][j4] += ej1[indx[i4]][jb] * tsmall[i4][j4];
            }
          libtetrabz_dblstep2(nb, ei2, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j4 = 0; j4 < 4; j4++)
                w1[ib][indx[i4]][jb] += v * tsmall[i4][j4] * w2[jb][j4];
        }

        libtetrabz_tsmall_b2(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (j4 = 0; j4 < 4; j4++) {
            ei2[j4] = 0.0;
            for (jb = 0; jb < nb; jb++) ej2[jb][j4] = 0.0;
          }
          for (i4 = 0; i4 < 4; i4++)
            for (j4 = 0; j4 < 4; j4++) {
              ei2[j4] += ei1[indx[i4]][ib] * tsmall[i4][j4];
              for (jb = 0; jb < nb; jb++)
                ej2[jb][j4] += ej1[indx[i4]][jb] * tsmall[i4][j4];
            }
          libtetrabz_dblstep2(nb, ei2, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j4 = 0; j4 < 4; j4++)
                w1[ib][indx[i4]][jb] += v * tsmall[i4][j4] * w2[jb][j4];
        }

        libtetrabz_tsmall_b3(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (j4 = 0; j4 < 4; j4++) {
            ei2[j4] = 0.0;
            for (jb = 0; jb < nb; jb++) ej2[jb][j4] = 0.0;
          }
          for (i4 = 0; i4 < 4; i4++)
            for (j4 = 0; j4 < 4; j4++) {
              ei2[j4] += ei1[indx[i4]][ib] * tsmall[i4][j4];
              for (jb = 0; jb < nb; jb++)
                ej2[jb][j4] += ej1[indx[i4]][jb] * tsmall[i4][j4];
            }
          libtetrabz_dblstep2(nb, ei2, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j4 = 0; j4 < 4; j4++)
                w1[ib][indx[i4]][jb] += v * tsmall[i4][j4] * w2[jb][j4];
        }
      }
      else if (e[2] <= 0.0 && 0.0 < e[3]) {

        libtetrabz_tsmall_c1(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (j4 = 0; j4 < 4; j4++) {
            ei2[j4] = 0.0;
            for (jb = 0; jb < nb; jb++) ej2[jb][j4] = 0.0;
          }
          for (i4 = 0; i4 < 4; i4++)
            for (j4 = 0; j4 < 4; j4++) {
              ei2[j4] += ei1[indx[i4]][ib] * tsmall[i4][j4];
              for (jb = 0; jb < nb; jb++)
                ej2[jb][j4] += ej1[indx[i4]][jb] * tsmall[i4][j4];
            }
          libtetrabz_dblstep2(nb, ei2, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j4 = 0; j4 < 4; j4++)
                w1[ib][indx[i4]][jb] += v * tsmall[i4][j4] * w2[jb][j4];
        }

        libtetrabz_tsmall_c2(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (j4 = 0; j4 < 4; j4++) {
            ei2[j4] = 0.0;
            for (jb = 0; jb < nb; jb++) ej2[jb][j4] = 0.0;
          }
          for (i4 = 0; i4 < 4; i4++)
            for (j4 = 0; j4 < 4; j4++) {
              ei2[j4] += ei1[indx[i4]][ib] * tsmall[i4][j4];
              for (jb = 0; jb < nb; jb++)
                ej2[jb][j4] += ej1[indx[i4]][jb] * tsmall[i4][j4];
            }
          libtetrabz_dblstep2(nb, ei2, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j4 = 0; j4 < 4; j4++)
                w1[ib][indx[i4]][jb] += v * tsmall[i4][j4] * w2[jb][j4];
        }

        libtetrabz_tsmall_c3(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (j4 = 0; j4 < 4; j4++) {
            ei2[j4] = 0.0;
            for (jb = 0; jb < nb; jb++) ej2[jb][j4] = 0.0;
          }
          for (i4 = 0; i4 < 4; i4++)
            for (j4 = 0; j4 < 4; j4++) {
              ei2[j4] += ei1[indx[i4]][ib] * tsmall[i4][j4];
              for (jb = 0; jb < nb; jb++)
                ej2[jb][j4] += ej1[indx[i4]][jb] * tsmall[i4][j4];
            }
          libtetrabz_dblstep2(nb, ei2, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j4 = 0; j4 < 4; j4++)
                w1[ib][indx[i4]][jb] += v * tsmall[i4][j4] * w2[jb][j4];
        }
      }
      else if (e[3] <= 0.0) {
        for (i4 = 0; i4 < 4; i4++) {
          ei2[i4] = ei1[i4][ib];
          for (jb = 0; jb < nb; jb++)
            ej2[jb][i4] = ej1[i4][jb];
        }
        libtetrabz_dblstep2(nb, ei2, ej2, w2);
        for (i4 = 0; i4 < 4; i4++)
          for (jb = 0; jb < nb; jb++)
            w1[ib][i4][jb] += w2[jb][i4];
      }
      else {
        continue;
      }
    }
    for (i20 = 0; i20 < 20; i20++)
      for (ib = 0; ib < nb; ib++)
        for (i4 = 0; i4 < 4; i4++)
          for (jb = 0; jb < nb; jb++)
            wght[ikv[it][i20]][ib][jb] += wlsm[i20][i4] * w1[ib][i4][jb];
  }
  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      for (jb = 0; jb < nb; jb++)
        wght[ik][ib][jb] /= (6.0 * (double) nk);

  free(ikv[0]);
  free(ikv);
  free(ei1[0]);
  free(ei1);
  free(ej1[0]);
  free(ej1);
  free(ej2[0]);
  free(ej2);
  free(w1[0][0]);
  free(w1[0]);
  free(w1);
  free(w2[0]);
  free(w2);
  /*
  Convert weight to python list object
  */
  wght_po = PyList_New(nk * nb * nb);
  ierr = 0;
  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      for (jb = 0; jb < nb; jb++)
        ierr = PyList_SetItem(wght_po, ik * nb * nb + ib * nb + jb, PyFloat_FromDouble(wght[ik][ib][jb]));
  if (ierr != 0)printf("Error in PyList_SetItem\n");

  free(eig1[0]);
  free(eig1);
  free(eig2[0]);
  free(eig2);
  free(wght[0][0]);
  free(wght[0]);
  free(wght);

  return wght_po;
}
/*
Compute Dos : Delta(E - E1)
*/
static PyObject* dos_c(PyObject* self, PyObject* args) {
  int it, ik, ib, ii, jj, ie, ** ikv, indx[4], ierr, ng[3], nk, nb, ne;
  double wlsm[20][4], ** ei1, e[4], *** w1, v, tsmall[4][3], bvec[3][3], ** eig, * e0, *** wght;
  PyObject* eig_po, * e0_po, * wght_po;
  /*
  Read input from python object
  */
  if (!PyArg_ParseTuple(args, "iiiiiidddddddddOO",
    &ng[0], &ng[1], &ng[2], &nk, &nb, &ne,
    &bvec[0][0], &bvec[0][1], &bvec[0][2], &bvec[1][0], &bvec[1][1], &bvec[1][2], &bvec[2][0], &bvec[2][1], &bvec[2][2],
    &eig_po, &e0_po))
    return NULL;
  /*
  convert python list object to array
  */
  eig = (double**)malloc(nk * sizeof(double*));
  eig[0] = (double*)malloc(nk * nb * sizeof(double));
  wght = (double***)malloc(nk * sizeof(double**));
  wght[0] = (double**)malloc(nk * nb * sizeof(double*));
  wght[0][0] = (double*)malloc(nk * nb * ne * sizeof(double));
  for (ik = 0; ik < nk; ik++) {
    eig[ik] = eig[0] + ik * nb;
    wght[ik] = wght[0] + ik * nb;
    for (ib = 0; ib < nb; ib++) {
      wght[ik][ib] = wght[0][0] + ik * nb * ne + ib * ne;
    }
  }
  e0 = (double*)malloc(ne * sizeof(double));

  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      eig[ik][ib] = PyFloat_AsDouble(PyList_GetItem(eig_po, ik * nb + ib));
  for (ie = 0; ie < ne; ie++)e0[ie] = PyFloat_AsDouble(PyList_GetItem(e0_po, ie));
  /*
  Start main calculation
  */
  ikv = (int**)malloc(6 * nk * sizeof(int*));
  ikv[0] = (int*)malloc(6 * nk * 20 * sizeof(int));
  for (ik = 0; ik < 6 * nk; ik++) {
    ikv[ik] = ikv[0] + ik * 20;
  }

  ei1 = (double**)malloc(4 * sizeof(double*));
  ei1[0] = (double*)malloc(4 * nb * sizeof(double));
  for (ii = 0; ii < 4; ii++) {
    ei1[ii] = ei1[0] + ii * nb;
  }

  w1 = (double***)malloc(nb * sizeof(double**));
  w1[0] = (double**)malloc(nb * ne * sizeof(double*));
  w1[0][0] = (double*)malloc(nb * ne * 4 * sizeof(double));
  for (ib = 0; ib < nb; ib++) {
    w1[ib] = w1[0] + ib * ne;
    for (ie = 0; ie < ne; ie++) {
      w1[ib][ie] = w1[0][0] + ib * ne * 4 + ie * 4;
    }
  }

  libtetrabz_initialize(ng, bvec, wlsm, ikv);

  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      for (ie = 0; ie < ne; ie++)
        wght[ik][ib][ie] = 0.0;

  for (it = 0; it < 6 * nk; it++) {

    for (ii = 0; ii < 4; ii++)
      for (ib = 0; ib < nb; ib++)
        ei1[ii][ib] = 0.0;
    for (jj = 0; jj < 20; jj++)
      for (ii = 0; ii < 4; ii++)
        for (ib = 0; ib < nb; ib++)
          ei1[ii][ib] += eig[ikv[it][jj]][ib] * wlsm[jj][ii];

    for (ib = 0; ib < nb; ib++) {

      for (ie = 0; ie < ne; ie++)
        for (ii = 0; ii < 4; ii++)
          w1[ib][ie][ii] = 0.0;

      for (ii = 0; ii < 4; ii++) e[ii] = ei1[ii][ib];
      eig_sort(4, e, indx);

      for (ie = 0; ie < ne; ie++) {

        if ((e[0] < e0[ie] && e0[ie] <= e[1]) || (e[0] <= e0[ie] && e0[ie] < e[1])) {

          libtetrabz_triangle_a1(e, e0[ie], &v, tsmall);
          for (ii = 0; ii < 4; ii++)
            for (jj = 0; jj < 3; jj++)
              w1[ib][ie][indx[ii]] += v * tsmall[ii][jj] / 3.0;

        }
        else if ((e[1] < e0[ie] && e0[ie] <= e[2]) || (e[1] <= e0[ie] && e0[ie] < e[2])) {

          libtetrabz_triangle_b1(e, e0[ie], &v, tsmall);
          for (ii = 0; ii < 4; ii++)
            for (jj = 0; jj < 3; jj++)
              w1[ib][ie][indx[ii]] += v * tsmall[ii][jj] / 3.0;

          libtetrabz_triangle_b2(e, e0[ie], &v, tsmall);
          for (ii = 0; ii < 4; ii++)
            for (jj = 0; jj < 3; jj++)
              w1[ib][ie][indx[ii]] += v * tsmall[ii][jj] / 3.0;
        }
        else if ((e[2] < e0[ie] && e0[ie] <= e[3]) || (e[2] <= e0[ie] && e0[ie] < e[3])) {

          libtetrabz_triangle_c1(e, e0[ie], &v, tsmall);
          for (ii = 0; ii < 4; ii++)
            for (jj = 0; jj < 3; jj++)
              w1[ib][ie][indx[ii]] += v * tsmall[ii][jj] / 3.0;
        }
        else {
          continue;
        }
      }
    }
    for (ii = 0; ii < 20; ii++)
      for (ib = 0; ib < nb; ib++)
        for (ie = 0; ie < ne; ie++)
          for (jj = 0; jj < 4; jj++)
            wght[ikv[it][ii]][ib][ie] += wlsm[ii][jj] * w1[ib][ie][jj];
  }
  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      for (ie = 0; ie < ne; ie++)
        wght[ik][ib][ie] /= (6.0 * (double)nk);

  free(ikv[0]);
  free(ikv);
  free(ei1[0]);
  free(ei1);
  free(w1[0][0]);
  free(w1[0]);
  free(w1);
  /*
  Convert weight to python list object
  */
  wght_po = PyList_New(nk * nb * ne);
  ierr = 0;
  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      for (ie = 0; ie < ne; ie++)
        ierr = PyList_SetItem(wght_po, ik * nb * ne + ib * ne + ie, PyFloat_FromDouble(wght[ik][ib][ie]));
  if (ierr != 0)printf("Error in PyList_SetItem\n");

  free(eig[0]);
  free(eig);
  free(wght[0][0]);
  free(wght[0]);
  free(wght);
  free(e0);

  return wght_po;
}
/*
Compute integrated Dos : theta(E - E1)
*/
static PyObject* intdos_c(PyObject* self, PyObject* args) {
  int it, ik, ib, ii, jj, ie, ** ikv, indx[4], ierr, ng[3], nk, nb, ne;
  double wlsm[20][4], ** ei1, e[4], *** w1, v, tsmall[4][4], bvec[3][3], ** eig, * e0, *** wght;
  PyObject* eig_po, * e0_po, * wght_po;
  /*
  Read input from python object
  */
  if (!PyArg_ParseTuple(args, "iiiiiidddddddddOO",
    &ng[0], &ng[1], &ng[2], &nk, &nb, &ne,
    &bvec[0][0], &bvec[0][1], &bvec[0][2], &bvec[1][0], &bvec[1][1], &bvec[1][2], &bvec[2][0], &bvec[2][1], &bvec[2][2],
    &eig_po, &e0_po))
    return NULL;
  /*
  convert python list object to array
  */
  eig = (double**)malloc(nk * sizeof(double*));
  eig[0] = (double*)malloc(nk * nb * sizeof(double));
  wght = (double***)malloc(nk * sizeof(double**));
  wght[0] = (double**)malloc(nk * nb * sizeof(double*));
  wght[0][0] = (double*)malloc(nk * nb * ne * sizeof(double));
  for (ik = 0; ik < nk; ik++) {
    eig[ik] = eig[0] + ik * nb;
    wght[ik] = wght[0] + ik * nb;
    for (ib = 0; ib < nb; ib++) {
      wght[ik][ib] = wght[0][0] + ik * nb * ne + ib * ne;
    }
  }
  e0 = (double*)malloc(ne * sizeof(double));

  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      eig[ik][ib] = PyFloat_AsDouble(PyList_GetItem(eig_po, ik * nb + ib));
  for (ie = 0; ie < ne; ie++)e0[ie] = PyFloat_AsDouble(PyList_GetItem(e0_po, ie));
  /*
  Start main calculation
  */
  ikv = (int**)malloc(6 * nk * sizeof(int*));
  ikv[0] = (int*)malloc(6 * nk * 20 * sizeof(int));
  for (ik = 0; ik < 6 * nk; ik++) {
    ikv[ik] = ikv[0] + ik * 20;
  }

  ei1 = (double**)malloc(4 * sizeof(double*));
  ei1[0] = (double*)malloc(4 * nb * sizeof(double));
  for (ii = 0; ii < 4; ii++) {
    ei1[ii] = ei1[0] + ii * nb;
  }

  w1 = (double***)malloc(nb * sizeof(double**));
  w1[0] = (double**)malloc(nb * ne * sizeof(double*));
  w1[0][0] = (double*)malloc(nb * ne * 4 * sizeof(double));
  for (ib = 0; ib < nb; ib++) {
    w1[ib] = w1[0] + ib * ne;
    for (ie = 0; ie < ne; ie++) {
      w1[ib][ie] = w1[0][0] + ib * ne * 4 + ie * 4;
    }
  }

  libtetrabz_initialize(ng, bvec, wlsm, ikv);

  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      for (ie = 0; ie < ne; ie++)
        wght[ik][ib][ie] = 0.0;

  for (it = 0; it < 6 * nk; it++) {

    for (ii = 0; ii < 4; ii++)
      for (ib = 0; ib < nb; ib++)
        ei1[ii][ib] = 0.0;
    for (jj = 0; jj < 20; jj++)
      for (ii = 0; ii < 4; ii++)
        for (ib = 0; ib < nb; ib++)
          ei1[ii][ib] += eig[ikv[it][jj]][ib] * wlsm[jj][ii];

    for (ib = 0; ib < nb; ib++)
      for (ie = 0; ie < ne; ie++)
        for (ii = 0; ii < 4; ii++)
          w1[ib][ie][ii] = 0.0;

    for (ib = 0; ib < nb; ib++) {

      for (ii = 0; ii < 4; ii++) e[ii] = ei1[ii][ib];
      eig_sort(4, e, indx);

      for (ie = 0; ie < ne; ie++) {

        if ((e[0] <= e0[ie] && e0[ie] < e[1]) || (e[0] < e0[ie] && e0[ie] <= e[1])) {
          libtetrabz_tsmall_a1(e, e0[ie], &v, tsmall);
          for (ii = 0; ii < 4; ii++)
            for (jj = 0; jj < 4; jj++)
              w1[ib][ie][indx[ii]] += v * tsmall[ii][jj] * 0.25;
        }
        else if ((e[1] <= e0[ie] && e0[ie] < e[2]) || (e[1] < e0[ie] && e0[ie] <= e[2])) {

          libtetrabz_tsmall_b1(e, e0[ie], &v, tsmall);
          for (ii = 0; ii < 4; ii++)
            for (jj = 0; jj < 4; jj++)
              w1[ib][ie][indx[ii]] += v * tsmall[ii][jj] * 0.25;

          libtetrabz_tsmall_b2(e, e0[ie], &v, tsmall);
          for (ii = 0; ii < 4; ii++)
            for (jj = 0; jj < 4; jj++)
              w1[ib][ie][indx[ii]] += v * tsmall[ii][jj] * 0.25;

          libtetrabz_tsmall_b3(e, e0[ie], &v, tsmall);
          for (ii = 0; ii < 4; ii++)
            for (jj = 0; jj < 4; jj++)
              w1[ib][ie][indx[ii]] += v * tsmall[ii][jj] * 0.25;

        }
        else if ((e[2] <= e0[ie] && e0[ie] < e[3]) || (e[2] < e0[ie] && e0[ie] <= e[3])) {

          libtetrabz_tsmall_c1(e, e0[ie], &v, tsmall);
          for (ii = 0; ii < 4; ii++)
            for (jj = 0; jj < 4; jj++)
              w1[ib][ie][indx[ii]] += v * tsmall[ii][jj] * 0.25;

          libtetrabz_tsmall_c2(e, e0[ie], &v, tsmall);
          for (ii = 0; ii < 4; ii++)
            for (jj = 0; jj < 4; jj++)
              w1[ib][ie][indx[ii]] += v * tsmall[ii][jj] * 0.25;

          libtetrabz_tsmall_c3(e, e0[ie], &v, tsmall);
          for (ii = 0; ii < 4; ii++)
            for (jj = 0; jj < 4; jj++)
              w1[ib][ie][indx[ii]] += v * tsmall[ii][jj] * 0.25;

        }
        else if (e[3] <= e0[ie]) {
          for (ii = 0; ii < 4; ii++)
            w1[ib][ie][ii] = 0.25;
        }
        else {
          continue;
        }
      }
    }

    for (ii = 0; ii < 20; ii++)
      for (ib = 0; ib < nb; ib++)
        for (ie = 0; ie < ne; ie++)
          for (jj = 0; jj < 4; jj++)
            wght[ikv[it][ii]][ib][ie] += wlsm[ii][jj] * w1[ib][ie][jj];
  }
  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      for (ie = 0; ie < ne; ie++)
        wght[ik][ib][ie] /= (6.0 * (double)nk);

  free(ikv[0]);
  free(ikv);
  free(ei1[0]);
  free(ei1);
  free(w1[0][0]);
  free(w1[0]);
  free(w1);
  /*
  Convert weight to python list object
  */
  wght_po = PyList_New(nk * nb * ne);
  ierr = 0;
  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      for (ie = 0; ie < ne; ie++)
        ierr = PyList_SetItem(wght_po, ik * nb * ne + ib * ne + ie, PyFloat_FromDouble(wght[ik][ib][ie]));
  if (ierr != 0)printf("Error in PyList_SetItem\n");

  free(eig[0]);
  free(eig);
  free(wght[0][0]);
  free(wght[0]);
  free(wght);
  free(e0);

  return wght_po; }
/*
 3rd step for Fermi's golden rule
*/
static void libtetrabz_fermigr3(
    int ne,
    double *e0,
    double de[4],
    double **w1
) {
  int i4, j3, ie, indx[4];
  double e[4], tsmall[4][3], v;

  for (i4 = 0; i4 < 4; i4++) e[i4] = de[i4];
  eig_sort(4, e, indx);

  for (ie = 0; ie < ne; ie++) {

    for (i4 = 0; i4 < 4; i4++) w1[i4][ie] = 0.0;

    if (e[0] < e0[ie] && e0[ie] <= e[1]) {

      libtetrabz_triangle_a1(e, e0[ie], &v, tsmall);
      for (i4 = 0; i4 < 4; i4++)
        for (j3 = 0; j3 < 3; j3++)
          w1[indx[i4]][ie] += v * tsmall[i4][j3] / 3.0;
    }
    else if (e[1] < e0[ie] && e0[ie] <= e[2]) {

      libtetrabz_triangle_b1(e, e0[ie], &v, tsmall);
      for (i4 = 0; i4 < 4; i4++)
        for (j3 = 0; j3 < 3; j3++)
          w1[indx[i4]][ie] += v * tsmall[i4][j3] / 3.0;

      libtetrabz_triangle_b2(e, e0[ie], &v, tsmall);
      for (i4 = 0; i4 < 4; i4++)
        for (j3 = 0; j3 < 3; j3++)
          w1[indx[i4]][ie] += v * tsmall[i4][j3] / 3.0;
    }
    else if (e[2] < e0[ie] && e0[ie] < e[3]) {

      libtetrabz_triangle_c1(e, e0[ie], &v, tsmall);
      for (i4 = 0; i4 < 4; i4++)
        for (j3 = 0; j3 < 3; j3++)
          w1[indx[i4]][ie] += v * tsmall[i4][j3] / 3.0;
    }
  }
}
/*
2nd step for Fermi's golden rule
*/
static void libtetrabz_fermigr2(
    int nb,
    int ne,
    double *e0,
    double *ei1,
    double **ej1,
    double ***w1
) {
  int ib, i4, j4, ie, indx[4];
  double e[4], tsmall[4][4], v, de[4], thr, **w2;

  w2 = (double**)malloc(4 * sizeof(double*));
  w2[0] = (double*)malloc(4 * ne * sizeof(double));
  for (i4 = 0; i4 < 4; i4++) {
    w2[i4] = w2[0] + i4 * ne;
  }

  thr = 1.0e-8;

  for (ib = 0; ib < nb; ib++) {

    for (i4 = 0; i4 < 4; i4++)
      for (ie = 0; ie < ne; ie++)
        w1[ib][i4][ie] = 0.0;

    for (i4 = 0; i4 < 4; i4++) e[i4] = -ej1[ib][i4];
    eig_sort(4, e, indx);

    if ((e[0] <= 0.0 && 0.0 < e[1]) || (e[0] < 0.0 && 0.0 <= e[1])) {

      libtetrabz_tsmall_a1(e, 0.0, &v, tsmall);

      if (v > thr) {
        for (j4 = 0; j4 < 4; j4++) de[j4] = 0.0;
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            de[j4] += (ej1[ib][indx[i4]] - ei1[indx[i4]]) * tsmall[i4][j4];
        libtetrabz_fermigr3(ne, e0, de, w2);
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            for (ie = 0; ie < ne; ie++)
              w1[ib][indx[i4]][ie] += v * tsmall[i4][j4] * w2[j4][ie];
      }
    }
    else if ((e[1] <= 0.0 && 0.0 < e[2]) || (e[1] < 0.0 && 0.0 <= e[2])) {

      libtetrabz_tsmall_b1(e, 0.0, &v, tsmall);

      if (v > thr) {
        for (j4 = 0; j4 < 4; j4++) de[j4] = 0.0;
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            de[j4] += (ej1[ib][indx[i4]] - ei1[indx[i4]]) * tsmall[i4][j4];
        libtetrabz_fermigr3(ne, e0, de, w2);
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            for (ie = 0; ie < ne; ie++)
              w1[ib][indx[i4]][ie] += v * tsmall[i4][j4] * w2[j4][ie];
      }

      libtetrabz_tsmall_b2(e, 0.0, &v, tsmall);

      if (v > thr) {
        for (j4 = 0; j4 < 4; j4++) de[j4] = 0.0;
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            de[j4] += (ej1[ib][indx[i4]] - ei1[indx[i4]]) * tsmall[i4][j4];
        libtetrabz_fermigr3(ne, e0, de, w2);
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            for (ie = 0; ie < ne; ie++)
              w1[ib][indx[i4]][ie] += v * tsmall[i4][j4] * w2[j4][ie];
      }

      libtetrabz_tsmall_b3(e, 0.0, &v, tsmall);

      if (v > thr) {
        for (j4 = 0; j4 < 4; j4++) de[j4] = 0.0;
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            de[j4] += (ej1[ib][indx[i4]] - ei1[indx[i4]]) * tsmall[i4][j4];
        libtetrabz_fermigr3(ne, e0, de, w2);
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            for (ie = 0; ie < ne; ie++)
              w1[ib][indx[i4]][ie] += v * tsmall[i4][j4] * w2[j4][ie];
      }
    }
    else if ((e[2] <= 0.0 && 0.0 < e[3]) || (e[2] < 0.0 && 0.0 <= e[3])) {

      libtetrabz_tsmall_c1(e, 0.0, &v, tsmall);
      
      if (v > thr) {
        for (j4 = 0; j4 < 4; j4++) de[j4] = 0.0;
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            de[j4] += (ej1[ib][indx[i4]] - ei1[indx[i4]]) * tsmall[i4][j4];
        libtetrabz_fermigr3(ne, e0, de, w2);
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            for (ie = 0; ie < ne; ie++)
              w1[ib][indx[i4]][ie] += v * tsmall[i4][j4] * w2[j4][ie];
      }

      libtetrabz_tsmall_c2(e, 0.0, &v, tsmall);

      if (v > thr) {
        for (j4 = 0; j4 < 4; j4++) de[j4] = 0.0;
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            de[j4] += (ej1[ib][indx[i4]] - ei1[indx[i4]]) * tsmall[i4][j4];
        libtetrabz_fermigr3(ne, e0, de, w2);
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            for (ie = 0; ie < ne; ie++)
              w1[ib][indx[i4]][ie] += v * tsmall[i4][j4] * w2[j4][ie];
      }

      libtetrabz_tsmall_c3(e, 0.0, &v, tsmall);

      if (v > thr) {
        for (j4 = 0; j4 < 4; j4++) de[j4] = 0.0;
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            de[j4] += (ej1[ib][indx[i4]] - ei1[indx[i4]]) * tsmall[i4][j4];
        libtetrabz_fermigr3(ne, e0, de, w2);
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            for (ie = 0; ie < ne; ie++)
              w1[ib][indx[i4]][ie] += v * tsmall[i4][j4] * w2[j4][ie];
      }
    }
    else if (e[3] <= 0.0) {
      for (i4 = 0; i4 < 4; i4++)
        de[i4] = ej1[ib][i4] - ei1[i4];
      libtetrabz_fermigr3(ne, e0, de, w2);
      for (i4 = 0; i4 < 4; i4++)
        for (ie = 0; ie < ne; ie++)
          w1[ib][i4][ie] += w2[i4][ie];
    }
  }

  free(w2[0]);
  free(w2);
}
/*
Main SUBROUTINE for Fermi's Golden rule : Theta(- E1) * Theta(E2) * Delta(E2 - E1 - w)
*/
static PyObject* fermigr_c(PyObject* self, PyObject* args)
{
  int it, ik, ib, i20, i4, j4, jb, ie, **ikv, indx[4], ierr, ng[3], nk, nb, ne;
  double wlsm[20][4], **ei1, **ej1, ei2[4], ** ej2, e[4], **** w1, *** w2, v, tsmall[4][4], thr,
    bvec[3][3], ** eig1, ** eig2, *e0, * ***wght;
  PyObject* eig1_po, * eig2_po, *e0_po, * wght_po;
  /*
  Read input from python object
  */
  if (!PyArg_ParseTuple(args, "iiiiiidddddddddOOO",
    &ng[0], &ng[1], &ng[2], &nk, &nb, &ne, 
    &bvec[0][0], &bvec[0][1], &bvec[0][2], &bvec[1][0], &bvec[1][1], &bvec[1][2], &bvec[2][0], &bvec[2][1], &bvec[2][2],
    &eig1_po, &eig2_po, &e0_po))
    return NULL;
  /*
  convert python list object to array
  */
  eig1 = (double**)malloc(nk * sizeof(double*));
  eig1[0] = (double*)malloc(nk * nb * sizeof(double));
  eig2 = (double**)malloc(nk * sizeof(double*));
  eig2[0] = (double*)malloc(nk * nb * sizeof(double));
  wght = (double****)malloc(nk * sizeof(double***));
  wght[0] = (double***)malloc(nk * nb * sizeof(double**));
  wght[0][0] = (double**)malloc(nk * nb * nb * sizeof(double*));
  wght[0][0][0] = (double*)malloc(nk * nb * nb *ne *sizeof(double));
  for (ik = 0; ik < nk; ik++) {
    eig1[ik] = eig1[0] + ik * nb;
    eig2[ik] = eig2[0] + ik * nb;
    wght[ik] = wght[0] + ik * nb;
    for (ib = 0; ib < nb; ib++) {
      wght[ik][ib] = wght[0][0] + ik * nb * nb + ib * nb;
      for (jb = 0; jb < nb; jb++) {
        wght[ik][ib][jb] = wght[0][0][0]
          + ik * nb * nb * ne + ib * nb * ne + jb * ne;
      }
    }
  }
  e0 = (double*)malloc(ne * sizeof(double));

  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++) {
      eig1[ik][ib] = PyFloat_AsDouble(PyList_GetItem(eig1_po, ik * nb + ib));
      eig2[ik][ib] = PyFloat_AsDouble(PyList_GetItem(eig2_po, ik * nb + ib));
    }
  for (ie = 0; ie < ne; ie++)e0[ie] = PyFloat_AsDouble(PyList_GetItem(e0_po, ie));
  /*
  Start main calculation
  */
  ikv = (int**)malloc(6 * nk * sizeof(int*));
  ikv[0] = (int*)malloc(6 * nk * 20 * sizeof(int));
  for (it = 0; it < 6 * nk; it++) {
    ikv[it] = ikv[0] + it * 20;
  }

  ei1 = (double**)malloc(4 * sizeof(double*));
  ej1 = (double**)malloc(4 * sizeof(double*));
  ei1[0] = (double*)malloc(4 * nb * sizeof(double));
  ej1[0] = (double*)malloc(4 * nb * sizeof(double));
  for (i4 = 0; i4 < 4; i4++) {
    ei1[i4] = ei1[0] + i4 * nb;
    ej1[i4] = ej1[0] + i4 * nb;
  }

  ej2 = (double**)malloc(nb * sizeof(double*));
  ej2[0] = (double*)malloc(nb * 4 * sizeof(double));
  for (ib = 0; ib < nb; ib++)
    ej2[ib] = ej2[0] + ib * 4;

  w1 = (double****)malloc(nb * sizeof(double***));
  w1[0] = (double***)malloc(nb * 4 * sizeof(double**));
  w1[0][0] = (double**)malloc(nb * 4 * nb * sizeof(double*));
  w1[0][0][0] = (double*)malloc(nb * 4 * nb * ne * sizeof(double));
  for (ib = 0; ib < nb; ib++) {
    w1[ib] = w1[0] + ib * 4;
    for (i4 = 0; i4 < 4; i4++) {
      w1[ib][i4] = w1[0][0] + ib * 4 * nb + i4 * nb;
      for (jb = 0; jb < nb; jb++) {
        w1[ib][i4][jb] = w1[0][0][0] + ib * 4 * nb * ne + i4 * nb * ne + jb * ne;
      }
    }
  }

  w2 = (double***)malloc(nb * sizeof(double**));
  w2[0] = (double**)malloc(nb * 4 * sizeof(double*));
  w2[0][0] = (double*)malloc(nb * 4 * ne * sizeof(double));
  for (ib = 0; ib < nb; ib++) {
    w2[ib] = w2[0] + ib * 4;
    for (i4 = 0; i4 < 4; i4++) {
      w2[ib][i4] = w2[0][0] + ib * 4 * ne + i4 * ne;
    }
  }

  libtetrabz_initialize(ng, bvec, wlsm, ikv);

  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      for (jb = 0; jb < nb; jb++)
        for (ie = 0; ie < ne; ie++)
          wght[ik][ib][jb][ie] = 0.0;

  thr = 1.0e-10;

  for (it = 0; it < 6*nk; it++) {

    for (i4 = 0; i4 < 4; i4++)
      for (ib = 0; ib < nb; ib++) {
        ei1[i4][ib] = 0.0;
        ej1[i4][ib] = 0.0;
      }
    for (i20 = 0; i20 < 20; i20++) {
      for (i4 = 0; i4 < 4; i4++) {
        for (ib = 0; ib < nb; ib++) {
          ei1[i4][ib] += eig1[ikv[it][i20]][ib] * wlsm[i20][i4];
          ej1[i4][ib] += eig2[ikv[it][i20]][ib] * wlsm[i20][i4];
        }
      }
    }
    
    for (ib = 0; ib < nb; ib++) {

      for (i4 = 0; i4 < 4; i4++)
        for (jb = 0; jb < nb; jb++)
          for (ie = 0; ie < ne; ie++)
            w1[ib][i4][jb][ie] = 0.0;

      for (i4 = 0; i4 < 4; i4++) e[i4] = ei1[i4][ib];
      eig_sort(4, e, indx);

      if (e[0] <= 0.0 && 0.0 < e[1]) {

        libtetrabz_tsmall_a1(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (j4 = 0; j4 < 4; j4++) {
            ei2[j4] = 0.0;
            for (jb = 0; jb < nb; jb++) ej2[jb][j4] = 0.0;
          }
          for (i4 = 0; i4 < 4; i4++)
            for (j4 = 0; j4 < 4; j4++) {
              ei2[j4] += ei1[indx[i4]][ib] * tsmall[i4][j4];
              for (jb = 0; jb < nb; jb++)
                ej2[jb][j4] += ej1[indx[i4]][jb] * tsmall[i4][j4];
            }
          libtetrabz_fermigr2(nb, ne, e0, ei2, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j4 = 0; j4 < 4; j4++)
                for (ie = 0; ie < ne; ie++)
                  w1[ib][indx[i4]][jb][ie] += v * tsmall[i4][j4] * w2[jb][j4][ie];
        }
      }
      else if (e[1] <= 0.0 && 0.0 < e[2]) {

        libtetrabz_tsmall_b1(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (j4 = 0; j4 < 4; j4++) {
            ei2[j4] = 0.0;
            for (jb = 0; jb < nb; jb++) ej2[jb][j4] = 0.0;
          }
          for (i4 = 0; i4 < 4; i4++)
            for (j4 = 0; j4 < 4; j4++) {
              ei2[j4] += ei1[indx[i4]][ib] * tsmall[i4][j4];
              for (jb = 0; jb < nb; jb++)
                ej2[jb][j4] += ej1[indx[i4]][jb] * tsmall[i4][j4];
            }
          libtetrabz_fermigr2(nb, ne, e0, ei2, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j4 = 0; j4 < 4; j4++)
                for (ie = 0; ie < ne; ie++)
                  w1[ib][indx[i4]][jb][ie] += v * tsmall[i4][j4] * w2[jb][j4][ie];
        }

        libtetrabz_tsmall_b2(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (j4 = 0; j4 < 4; j4++) {
            ei2[j4] = 0.0;
            for (jb = 0; jb < nb; jb++) ej2[jb][j4] = 0.0;
          }
          for (i4 = 0; i4 < 4; i4++)
            for (j4 = 0; j4 < 4; j4++) {
              ei2[j4] += ei1[indx[i4]][ib] * tsmall[i4][j4];
              for (jb = 0; jb < nb; jb++)
                ej2[jb][j4] += ej1[indx[i4]][jb] * tsmall[i4][j4];
            }
          libtetrabz_fermigr2(nb, ne, e0, ei2, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j4 = 0; j4 < 4; j4++)
                for (ie = 0; ie < ne; ie++)
                  w1[ib][indx[i4]][jb][ie] += v * tsmall[i4][j4] * w2[jb][j4][ie];
        }

        libtetrabz_tsmall_b3(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (j4 = 0; j4 < 4; j4++) {
            ei2[j4] = 0.0;
            for (jb = 0; jb < nb; jb++) ej2[jb][j4] = 0.0;
          }
          for (i4 = 0; i4 < 4; i4++)
            for (j4 = 0; j4 < 4; j4++) {
              ei2[j4] += ei1[indx[i4]][ib] * tsmall[i4][j4];
              for (jb = 0; jb < nb; jb++)
                ej2[jb][j4] += ej1[indx[i4]][jb] * tsmall[i4][j4];
            }
          libtetrabz_fermigr2(nb, ne, e0, ei2, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j4 = 0; j4 < 4; j4++)
                for (ie = 0; ie < ne; ie++)
                  w1[ib][indx[i4]][jb][ie] += v * tsmall[i4][j4] * w2[jb][j4][ie];
        }
      }
      else if (e[2] <= 0.0 && 0.0 < e[3]) {

        libtetrabz_tsmall_c1(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (j4 = 0; j4 < 4; j4++) {
            ei2[j4] = 0.0;
            for (jb = 0; jb < nb; jb++) ej2[jb][j4] = 0.0;
          }
          for (i4 = 0; i4 < 4; i4++)
            for (j4 = 0; j4 < 4; j4++) {
              ei2[j4] += ei1[indx[i4]][ib] * tsmall[i4][j4];
              for (jb = 0; jb < nb; jb++)
                ej2[jb][j4] += ej1[indx[i4]][jb] * tsmall[i4][j4];
            }
          libtetrabz_fermigr2(nb, ne, e0, ei2, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j4 = 0; j4 < 4; j4++)
                for (ie = 0; ie < ne; ie++)
                  w1[ib][indx[i4]][jb][ie] += v * tsmall[i4][j4] * w2[jb][j4][ie];
        }

        libtetrabz_tsmall_c2(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (j4 = 0; j4 < 4; j4++) {
            ei2[j4] = 0.0;
            for (jb = 0; jb < nb; jb++) ej2[jb][j4] = 0.0;
          }
          for (i4 = 0; i4 < 4; i4++)
            for (j4 = 0; j4 < 4; j4++) {
              ei2[j4] += ei1[indx[i4]][ib] * tsmall[i4][j4];
              for (jb = 0; jb < nb; jb++)
                ej2[jb][j4] += ej1[indx[i4]][jb] * tsmall[i4][j4];
            }
          libtetrabz_fermigr2(nb, ne, e0, ei2, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j4 = 0; j4 < 4; j4++)
                for (ie = 0; ie < ne; ie++)
                  w1[ib][indx[i4]][jb][ie] += v * tsmall[i4][j4] * w2[jb][j4][ie];
        }

        libtetrabz_tsmall_c3(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (j4 = 0; j4 < 4; j4++) {
            ei2[j4] = 0.0;
            for (jb = 0; jb < nb; jb++) ej2[jb][j4] = 0.0;
          }
          for (i4 = 0; i4 < 4; i4++)
            for (j4 = 0; j4 < 4; j4++) {
              ei2[j4] += ei1[indx[i4]][ib] * tsmall[i4][j4];
              for (jb = 0; jb < nb; jb++)
                ej2[jb][j4] += ej1[indx[i4]][jb] * tsmall[i4][j4];
            }
          libtetrabz_fermigr2(nb, ne, e0, ei2, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j4 = 0; j4 < 4; j4++)
                for (ie = 0; ie < ne; ie++)
                  w1[ib][indx[i4]][jb][ie] += v * tsmall[i4][j4] * w2[jb][j4][ie];
        }
      }
      else if (e[3] <= 0.0) {
        for (i4 = 0; i4 < 4; i4++) {
          ei2[i4] = ei1[i4][ib];
          for (jb = 0; jb < nb; jb++)
            ej2[jb][i4] = ej1[i4][jb];
        }
        libtetrabz_fermigr2(nb, ne, e0, ei2, ej2, w2);
        for (i4 = 0; i4 < 4; i4++)
          for (jb = 0; jb < nb; jb++)
            for (ie = 0; ie < ne; ie++)
              w1[ib][i4][jb][ie] += w2[jb][i4][ie];
      }
      else {
        continue;
      }
    }
    for (i20 = 0; i20 < 20; i20++)
      for (ib = 0; ib < nb; ib++)
        for (i4 = 0; i4 < 4; i4++)
          for (jb = 0; jb < nb; jb++)
            for (ie = 0; ie < ne; ie++)
              wght[ikv[it][i20]][ib][jb][ie] += wlsm[i20][i4] * w1[ib][i4][jb][ie];
  }
  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      for (jb = 0; jb < nb; jb++)
        for (ie = 0; ie < ne; ie++)
          wght[ik][ib][jb][ie] /= (6.0 * (double) nk);

  free(ikv[0]);
  free(ikv);
  free(ei1[0]);
  free(ei1);
  free(ej1[0]);
  free(ej1);
  free(ej2[0]);
  free(ej2);
  free(w1[0][0][0]);
  free(w1[0][0]);
  free(w1[0]);
  free(w1);
  free(w2[0][0]);
  free(w2[0]);
  free(w2);
  /*
  Convert weight to python list object
  */
  wght_po = PyList_New(nk * nb * nb * ne);
  ierr = 0;
  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      for (jb = 0; jb < nb; jb++)
        for (ie = 0; ie < ne; ie++)
          ierr = PyList_SetItem(wght_po, ik * nb * nb * ne + ib * nb * ne + jb * ne + ie,
            PyFloat_FromDouble(wght[ik][ib][jb][ie]));
  if (ierr != 0)printf("Error in PyList_SetItem\n");

  free(eig1[0]);
  free(eig1);
  free(eig2[0]);
  free(eig2);
  free(wght[0][0][0]);
  free(wght[0][0]);
  free(wght[0]);
  free(wght);
  free(e0);

  return wght_po;
}
/*
Main SUBROUTINE for occupation : Theta(EF - E1)
*/
void occ_main(
  int nk,
  int nb,
  double** eig,
  int **ikv,
  double wlsm[20][4],
  double **ei1,
  double **w1,
  double **wght
) {
  int it, ik, ib, ii, jj, indx[4];
  double e[4], v, tsmall[4][4];

  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      wght[ik][ib] = 0.0;

  for (it = 0; it < 6 * nk; it++) {

    for (ii = 0; ii < 4; ii++)
      for (ib = 0; ib < nb; ib++) 
        ei1[ii][ib] = 0.0;
    for (jj = 0; jj < 20; jj++) 
      for (ii = 0; ii < 4; ii++) 
        for (ib = 0; ib < nb; ib++) 
          ei1[ii][ib] += eig[ikv[it][jj]][ib] * wlsm[jj][ii];

    for (ib = 0; ib < nb; ib++) {

      for (ii = 0; ii < 4; ii++)
        w1[ib][ii] = 0.0;

      for (ii = 0; ii < 4; ii++) e[ii] = ei1[ii][ib];
      eig_sort(4, e, indx);

      if (e[0] <= 0.0 && 0.0 < e[1]) {
        libtetrabz_tsmall_a1(e, 0.0, &v, tsmall);
        for (ii = 0; ii < 4; ii++)
          for (jj = 0; jj < 4; jj++)
            w1[ib][indx[ii]] += v * tsmall[ii][jj] * 0.25;
      }
      else if (e[1] <= 0.0 && 0.0 < e[2]) {

        libtetrabz_tsmall_b1(e, 0.0, &v, tsmall);
        for (ii = 0; ii < 4; ii++)
          for (jj = 0; jj < 4; jj++)
            w1[ib][indx[ii]] += v * tsmall[ii][jj] * 0.25;

        libtetrabz_tsmall_b2(e, 0.0, &v, tsmall);
        for (ii = 0; ii < 4; ii++)
          for (jj = 0; jj < 4; jj++)
            w1[ib][indx[ii]] += v * tsmall[ii][jj] * 0.25;

        libtetrabz_tsmall_b3(e, 0.0, &v, tsmall);
        for (ii = 0; ii < 4; ii++)
          for (jj = 0; jj < 4; jj++)
            w1[ib][indx[ii]] += v * tsmall[ii][jj] * 0.25;
      }
      else if (e[2] <= 0.0 && 0.0 < e[3]) {

        libtetrabz_tsmall_c1(e, 0.0, &v, tsmall);
        for (ii = 0; ii < 4; ii++)
          for (jj = 0; jj < 4; jj++)
            w1[ib][indx[ii]] += v * tsmall[ii][jj] * 0.25;

        libtetrabz_tsmall_c2(e, 0.0, &v, tsmall);
        for (ii = 0; ii < 4; ii++)
          for (jj = 0; jj < 4; jj++)
            w1[ib][indx[ii]] += v * tsmall[ii][jj] * 0.25;

        libtetrabz_tsmall_c3(e, 0.0, &v, tsmall);
        for (ii = 0; ii < 4; ii++)
          for (jj = 0; jj < 4; jj++)
            w1[ib][indx[ii]] += v * tsmall[ii][jj] * 0.25;
      }
      else if (e[3] <= 0.0) {
        for (ii = 0; ii < 4; ii++)
          w1[ib][ii] += 0.25;
      }
      else {
        continue;
      }
    }
    for (ii = 0; ii < 20; ii++)
      for (ib = 0; ib < nb; ib++)
        for (jj = 0; jj < 4; jj++) {
          wght[ikv[it][ii]][ib] += wlsm[ii][jj] * w1[ib][jj];
        }
  }
  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      wght[ik][ib] /= (6.0 * (double)nk);
}
/*

*/
static PyObject* occ_c(PyObject* self, PyObject* args)
{
  int ik, ib, ii, nk, nb, ng[3], ierr, ** ikv;
  double** eig, ** wght, bvec[3][3], wlsm[20][4], **ei1, **w1;
  PyObject* eig_po, * wght_po;
  /*
  Read input from python object
  */
  if (!PyArg_ParseTuple(args, "iiiiidddddddddO",
    &ng[0], &ng[1], &ng[2], &nk, &nb,
    &bvec[0][0], &bvec[0][1], &bvec[0][2], &bvec[1][0], &bvec[1][1], &bvec[1][2], &bvec[2][0], &bvec[2][1], &bvec[2][2],
    &eig_po))
    return NULL;
  /*
  convert python list object to array
  */
  eig = (double**)malloc(nk * sizeof(double*));
  eig[0] = (double*)malloc(nk * nb * sizeof(double));
  wght = (double**)malloc(nk * sizeof(double*));
  wght[0] = (double*)malloc(nk * nb * sizeof(double));
  for (ik = 0; ik < nk; ik++) {
    eig[ik] = eig[0] + ik * nb;
    wght[ik] = wght[0] + ik * nb;
  }

  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      eig[ik][ib] = PyFloat_AsDouble(PyList_GetItem(eig_po, ik * nb + ib));
  
  /*
  Start main calculation
  */
  ikv = (int**)malloc(6 * nk * sizeof(int*));
  ikv[0] = (int*)malloc(6 * nk * 20 * sizeof(int));
  for (ik = 0; ik < 6 * nk; ik++) {
    ikv[ik] = ikv[0] + ik * 20;
  }
  ei1 = (double**)malloc(4 * sizeof(double*));
  ei1[0] = (double*)malloc(4 * nb * sizeof(double));
  for (ii = 0; ii < 4; ii++) {
    ei1[ii] = ei1[0] + ii * nb;
  }

  w1 = (double**)malloc(nb * sizeof(double*));
  w1[0] = (double*)malloc(nb * 4 * sizeof(double));
  for (ib = 0; ib < nb; ib++) {
    w1[ib] = w1[0] + ib * 4;
  }

  libtetrabz_initialize(ng, bvec, wlsm, ikv);

  occ_main(nk, nb, eig, ikv, wlsm, ei1, w1, wght);

  free(ikv[0]);
  free(ikv);
  free(ei1[0]);
  free(ei1);
  free(w1[0]);
  free(w1);
  /*
  Convert weight to python list object
  */
  wght_po = PyList_New(nk * nb);
  ierr = 0;
  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      ierr = PyList_SetItem(wght_po, ik * nb + ib, PyFloat_FromDouble(wght[ik][ib]));
  if (ierr != 0)printf("Error in PyList_SetItem\n");
  
  free(eig[0]);
  free(eig);
  free(wght[0]);
  free(wght);
  return wght_po;
}
/*
Calculate Fermi energy
*/
static PyObject* fermieng_c(PyObject* self, PyObject* args)
{
  int maxiter, ik, ib, iteration, nk, nb, ng[3],
    **ikv, ii, ierr;
  double eps, elw, eup, **eig2, sumkmid, bvec[3][3],
    **eig, **wght, wlsm[20][4], **ei1, **w1, nelec, ef;
  PyObject* eig_po, * wght_po;
  /*
  Read input from python object
  */
  if (!PyArg_ParseTuple(args, "iiiiidddddddddOd",
    &ng[0], &ng[1], &ng[2], &nk, &nb,
    &bvec[0][0], &bvec[0][1], &bvec[0][2], &bvec[1][0], &bvec[1][1], &bvec[1][2], &bvec[2][0], &bvec[2][1], &bvec[2][2],
    &eig_po, &nelec))
    return NULL;
  /*
  convert python list object to array
  */
  eig = (double**)malloc(nk * sizeof(double*));
  eig[0] = (double*)malloc(nk * nb * sizeof(double));
  wght = (double**)malloc(nk * sizeof(double*));
  wght[0] = (double*)malloc(nk * nb * sizeof(double));
  for (ik = 0; ik < nk; ik++) {
    eig[ik] = eig[0] + ik * nb;
    wght[ik] = wght[0] + ik * nb;
  }

  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      eig[ik][ib] = PyFloat_AsDouble(PyList_GetItem(eig_po, ik * nb + ib));

  /*
  Start main calculation
  */
  ikv = (int**)malloc(6 * nk * sizeof(int*));
  ikv[0] = (int*)malloc(6 * nk * 20 * sizeof(int));
  for (ik = 0; ik < 6 * nk; ik++) {
    ikv[ik] = ikv[0] + ik * 20;
  }
  ei1 = (double**)malloc(4 * sizeof(double*));
  ei1[0] = (double*)malloc(4 * nb * sizeof(double));
  for (ii = 0; ii < 4; ii++) {
    ei1[ii] = ei1[0] + ii * nb;
  }

  w1 = (double**)malloc(nb * sizeof(double*));
  w1[0] = (double*)malloc(nb * 4 * sizeof(double));
  for (ib = 0; ib < nb; ib++) {
    w1[ib] = w1[0] + ib * 4;
  }

  libtetrabz_initialize(ng, bvec, wlsm, ikv);

  eig2 = (double**)malloc(nk * sizeof(double*));
  eig2[0] = (double*)malloc(nk * nb * sizeof(double));
  for (ik = 0; ik < nk; ik++) {
    eig2[ik] = eig2[0] + ik * nb;
  }

  maxiter = 300;
  eps = 1.0e-10;

  elw = eig[0][0];
  eup = eig[0][0];
  for (ik = 0; ik < nk; ik++) {
    for (ib = 0; ib < nb; ib++) {
      if (elw > eig[ik][ib]) elw = eig[ik][ib];
      if (eup < eig[ik][ib]) eup = eig[ik][ib];
    }
  }
  /*
  Bisection method
  */
  for (iteration = 0; iteration < maxiter; iteration++) {

    ef = (eup + elw) * 0.5;
    /*
    Calc. # of electrons
    */
    for (ik = 0; ik < nk; ik++)
      for (ib = 0; ib < nb; ib++)
        eig2[ik][ib] = eig[ik][ib] - ef;
    occ_main(nk, nb, eig2, ikv, wlsm, ei1, w1, wght);

    sumkmid = 0.0;
    for (ik = 0; ik < nk; ik++) {
      for (ib = 0; ib < nb; ib++) {
        sumkmid += wght[ik][ib];
      }
    }
    /*
    convergence check
    */
    if (fabs(sumkmid - nelec) < eps) {
      break;
    }
    else if (sumkmid < nelec)
      elw = ef;
    else
      eup = ef;
  }
  free(eig2);
  free(ikv[0]);
  free(ikv);
  free(ei1[0]);
  free(ei1);
  free(w1[0]);
  free(w1);
  /*
  Convert weight to python list object
  */
  wght_po = PyList_New(nk * nb+2);
  ierr = 0;
  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      ierr = PyList_SetItem(wght_po, ik * nb + ib, PyFloat_FromDouble(wght[ik][ib]));
  ierr = PyList_SetItem(wght_po, nk*nb, PyFloat_FromDouble(ef));
  ierr = PyList_SetItem(wght_po, nk * nb+1, PyLong_FromLong((long)iteration));
  if (ierr != 0)printf("Error in PyList_SetItem\n");

  free(eig[0]);
  free(eig);
  free(wght[0]);
  free(wght);
  return wght_po;
}
/*
Results of Integration (1-x-y-z)/(g0+(g1-g0)x+(g2-g0)y+(g3-g0))
for 0<x<1, 0<y<1-x, 0<z<1-x-y
*/
/*
1, Different each other
*/
static void libtetrabz_polcmplx_1234(
  double g1,
  double g2,
  double g3,
  double g4,
  double* w
) {
  double w2, w3, w4;
  /*
  Real
  */
  w2 = 2.0 * (3.0 * g2 * g2 - 1.0) * (atan(g2) - atan(g1)) +
    (g2 * g2 - 3.0) * g2 * log((1.0 + g2 * g2) / (1.0 + g1 * g1));
  w2 = -2.0 * (g2 * g2 - 1.0) + w2 / (g2 - g1);
  w2 = w2 / (g2 - g1);
  w3 = 2.0 * (3.0 * g3 * g3 - 1.0) * (atan(g3) - atan(g1)) +
    (g3 * g3 - 3.0) * g3 * log((1.0 + g3 * g3) / (1.0 + g1 * g1));
  w3 = -2.0 * (g3 * g3 - 1.0) + w3 / (g3 - g1);
  w3 = w3 / (g3 - g1);
  w4 = 2.0 * (3.0 * g4 * g4 - 1.0) * (atan(g4) - atan(g1)) +
    (g4 * g4 - 3.0) * g4 * log((1.0 + g4 * g4) / (1.0 + g1 * g1));
  w4 = -2.0 * (g4 * g4 - 1.0) + w4 / (g4 - g1);
  w4 = w4 / (g4 - g1);
  w2 = (w2 - w3) / (g2 - g3);
  w4 = (w4 - w3) / (g4 - g3);
  w[0] = (w4 - w2) / (2.0 * (g4 - g2));
  /*
  Imaginal
  */
  w2 = 2.0 * (3.0 - g2 * g2) * g2 * (atan(g2) - atan(g1)) +
    (3.0 * g2 * g2 - 1.0) * log((1.0 + g2 * g2) / (1.0 + g1 * g1));
  w2 = 4.0 * g2 - w2 / (g2 - g1);
  w2 = w2 / (g2 - g1);
  w3 = 2.0 * (3.0 - g3 * g3) * g3 * (atan(g3) - atan(g1)) +
    (3.0 * g3 * g3 - 1.0) * log((1.0 + g3 * g3) / (1.0 + g1 * g1));
  w3 = 4.0 * g3 - w3 / (g3 - g1);
  w3 = w3 / (g3 - g1);
  w4 = 2.0 * (3.0 - g4 * g4) * g4 * (atan(g4) - atan(g1)) +
    (3.0 * g4 * g4 - 1.0) * log((1.0 + g4 * g4) / (1.0 + g1 * g1));
  w4 = 4.0 * g4 - w4 / (g4 - g1);
  w4 = w4 / (g4 - g1);
  w2 = (w2 - w3) / (g2 - g3);
  w4 = (w4 - w3) / (g4 - g3);
  w[1] = (w4 - w2) / (2.0 * (g4 - g2));
}
/*
2, g4 = g1
*/
static void libtetrabz_polcmplx_1231(
  double g1,
  double g2,
  double g3,
  double* w
) {
  double w2, w3;
  /*
  Real
  */
  w2 = 2.0 * (-1.0 + 3.0 * g2 * g2) * (atan(g2) - atan(g1))
    + g2 * (-3.0 + g2 * g2) * log((1.0 + g2 * g2) / (1.0 + g1 * g1));
  w2 = 2.0 * (1.0 - g2 * g2) + w2 / (g2 - g1);
  w2 = -g1 + w2 / (g2 - g1);
  w2 = w2 / (g2 - g1);
  w3 = 2.0 * (-1.0 + 3.0 * g3 * g3) * (atan(g3) - atan(g1))
    + g3 * (-3.0 + g3 * g3) * log((1.0 + g3 * g3) / (1.0 + g1 * g1));
  w3 = 2.0 * (1 - g3 * g3) + w3 / (g3 - g1);
  w3 = -g1 + w3 / (g3 - g1);
  w3 = w3 / (g3 - g1);
  w[0] = (w3 - w2) / (2.0 * (g3 - g2));
  /*
  Imaginal
  */
  w2 = 2.0 * g2 * (3.0 - g2 * g2) * (atan(g2) - atan(g1)) +
    (-1.0 + 3.0 * g2 * g2) * log((1.0 + g2 * g2) / (1.0 + g1 * g1));
  w2 = 4.0 * g2 - w2 / (g2 - g1);
  w2 = 1 + w2 / (g2 - g1);
  w2 = w2 / (g2 - g1);
  w3 = 2.0 * g3 * (3.0 - g3 * g3) * (atan(g3) - atan(g1)) +
    (-1.0 + 3.0 * g3 * g3) * log((1.0 + g3 * g3) / (1.0 + g1 * g1));
  w3 = 4.0 * g3 - w3 / (g3 - g1);
  w3 = 1 + w3 / (g3 - g1);
  w3 = w3 / (g3 - g1);
  w[1] = (w3 - w2) / (2.0 * (g3 - g2));
}
/*
3, g4 = g3
*/
static void libtetrabz_polcmplx_1233(
  double g1,
  double g2,
  double g3,
  double* w
) {
  double w2, w3;
  /*
  Real
  */
  w2 = 2.0 * (1.0 - 3.0 * g2 * g2) * (atan(g2) - atan(g1)) +
    g2 * (3.0 - g2 * g2) * log((1.0 + g2 * g2) / (1.0 + g1 * g1));
  w2 = 2.0 * (1 - g2 * g2) - w2 / (g2 - g1);
  w2 = w2 / (g2 - g1);
  w3 = 2.0 * (1.0 - 3.0 * g3 * g3) * (atan(g3) - atan(g1)) +
    g3 * (3.0 - g3 * g3) * log((1.0 + g3 * g3) / (1.0 + g1 * g1));
  w3 = 2.0 * (1 - g3 * g3) - w3 / (g3 - g1);
  w3 = w3 / (g3 - g1);
  w2 = (w3 - w2) / (g3 - g2);
  w3 = 4.0 * (1.0 - 3.0 * g1 * g3) * (atan(g3) - atan(g1))
    + (3.0 * g1 + 3.0 * g3 - 3.0 * g1 * g3 * g3 + g3 * g3 * g3) * log((1.0 + g3 * g3) / (1.0 + g1 * g1));
  w3 = -4.0 * (1.0 - g1 * g1) + w3 / (g3 - g1);
  w3 = 4.0 * g1 + w3 / (g3 - g1);
  w3 = w3 / (g3 - g1);
  w[0] = (w3 - w2) / (2.0 * (g3 - g2));
  /*
  Imaginal
  */
  w2 = 2.0 * g2 * (3.0 - g2 * g2) * (atan(g2) - atan(g1)) +
    (-1.0 + 3.0 * g2 * g2) * log((1.0 + g2 * g2) / (1.0 + g1 * g1));
  w2 = 4.0 * g2 - w2 / (g2 - g1);
  w2 = w2 / (g2 - g1);
  w3 = 2.0 * g3 * (3.0 - g3 * g3) * (atan(g3) - atan(g1)) +
    (-1.0 + 3.0 * g3 * g3) * log((1.0 + g3 * g3) / (1.0 + g1 * g1));
  w3 = 4.0 * g3 - w3 / (g3 - g1);
  w3 = w3 / (g3 - g1);
  w2 = (w3 - w2) / (g3 - g2);
  w3 = (3.0 * g1 - 3.0 * g1 * g3 * g3 + 3.0 * g3 + g3 * g3 * g3) * (atan(g3) - atan(g1))
    + (3.0 * g1 * g3 - 1.0) * log((1.0 + g3 * g3) / (1.0 + g1 * g1));
  w3 = w3 / (g3 - g1) - 4.0 * g1;
  w3 = w3 / (g3 - g1) - 2.0;
  w3 = (2.0 * w3) / (g3 - g1);
  w[1] = (w3 - w2) / (2.0 * (g3 - g2));
}
/*
4, g4 = g1 and g3 = g2
*/
static void libtetrabz_polcmplx_1221(
  double g1,
  double g2,
  double* w
) {
  /*
  Real
  */
  w[0] = -2.0 * (-1.0 + 2.0 * g1 * g2 + g2 * g2) * (atan(g2) - atan(g1))
    + (g1 + 2.0 * g2 - g1 * g2 * g2) * log((1.0 + g2 * g2) / (1.0 + g1 * g1));
  w[0] = 2.0 * (-1.0 + g1 * g1) + w[0] / (g2 - g1);
  w[0] = 3.0 * g1 + w[0] / (g2 - g1);
  w[0] = 2.0 + (3.0 * w[0]) / (g2 - g1);
  w[0] = w[0] / (2.0 * (g2 - g1));
  /*
  Imaginal
  */
  w[1] = 2.0 * (g1 + 2.0 * g2 - g1 * g2 * g2) * (atan(g2) - atan(g1))
    + (-1.0 + 2.0 * g1 * g2 + g2 * g2) * log((1 + g2 * g2) / (1 + g1 * g1));
  w[1] = -4.0 * g1 + w[1] / (g2 - g1);
  w[1] = -3.0 + w[1] / (g2 - g1);
  w[1] = (3.0 * w[1]) / (2.0 * (g2 - g1) * (g2 - g1));
}
/*
5, g4 = g3 = g2
*/
static void libtetrabz_polcmplx_1222(
  double g1,
  double g2,
  double* w
) {
  /*
  Real
  */
  w[0] = 2.0 * (-1.0 + g1 * g1 + 2.0 * g1 * g2) * (atan(g2) - atan(g1))
    + (-2.0 * g1 - g2 + g1 * g1 * g2) * log((1.0 + g2 * g2) / (1.0 + g1 * g1));
  w[0] = 2.0 * (1.0 - g1 * g1) + w[0] / (g2 - g1);
  w[0] = g1 - w[0] / (g2 - g1);
  w[0] = 1.0 - (3.0 * w[0]) / (g2 - g1);
  w[0] = w[0] / (2.0 * (g2 - g1));
  /*
  Imaginal
  */
  w[1] = 2.0 * (-2.0 * g1 - g2 + g1 * g1 * g2) * (atan(g2) - atan(g1))
    + (1.0 - g1 * g1 - 2.0 * g1 * g2) * log((1.0 + g2 * g2) / (1.0 + g1 * g1));
  w[1] = 4.0 * g1 + w[1] / (g2 - g1);
  w[1] = 1.0 + w[1] / (g2 - g1);
  w[1] = (3.0 * w[1]) / (2.0 * (g2 - g1) * (g2 - g1));
}
/*
6, g4 = g3 = g1
*/
static void libtetrabz_polcmplx_1211(
  double g1,
  double g2,
  double* w
) {
  /*
  Real
  */
  w[0] = 2.0 * (3.0 * g2 * g2 - 1.0) * (atan(g2) - atan(g1))
    + g2 * (g2 * g2 - 3.0) * log((1.0 + g2 * g2) / (1.0 + g1 * g1));
  w[0] = 2.0 * (1.0 - g1 * g1) + w[0] / (g2 - g1);
  w[0] = -5.0 * g1 + w[0] / (g2 - g1);
  w[0] = -11.0 + (3.0 * w[0]) / (g2 - g1);
  w[0] = w[0] / (6.0 * (g2 - g1));
  /*
  Imaginal
  */
  w[1] = 2.0 * g2 * (-3.0 + g2 * g2) * (atan(g2) - atan(g1))
    + (1.0 - 3.0 * g2 * g2) * log((1.0 + g2 * g2) / (1.0 + g1 * g1));
  w[1] = 4.0 * g2 + w[1] / (g2 - g1);
  w[1] = 1.0 + w[1] / (g2 - g1);
  w[1] = w[1] / (2.0 * (g2 - g1) * (g2 - g1));
}
/*
Tetrahedron method for delta(om - ep + e)
*/
static void libtetrabz_polcmplx3(
  int ne,
  double** e0,
  double de[4],
  double*** w1
) {
  int i4, ir, indx[4], ie;
  double e[4], x[4], thr, w2[4][2], denom;

  for (i4 = 0; i4 < 3; i4++) e[i4] = de[i4];
  eig_sort(4, e, indx);

  for (ie = 0; ie < ne; ie++) {
    /*
    I am not sure which one is better.
    The former is more stable.
    The latter is more accurate ?
    */
    for (i4 = 0; i4 < 4; i4++) {
      denom = (de[i4] + e0[ie][0]) * (de[i4] + e0[ie][0]) + e0[ie][1] * e0[ie][1];
      w1[i4][ie][0] = 0.25 * (de[i4] + e0[ie][0]) / denom;
      w1[i4][ie][1] = 0.25 * (-e0[ie][1]) / denom;
    }
    continue;

    for (i4 = 0; i4 < 4; i4++)
      x[i4] = (e[i4] + e0[ie][0]) / e0[ie][1];
    /* thr = maxval(de(1:4)) * 1d-3*/
    thr = fabs(x[0]);
    for (i4 = 0; i4 < 4; i4++)
      if (thr < fabs(x[i4])) thr = fabs(x[i4]);
    thr = fmax(1.0e-3, thr * 1.0e-2);

    if (fabs(x[3] - x[2]) < thr) {
      if (fabs(x[3] - x[1]) < thr) {
        if (fabs(x[3] - x[0]) < thr) {
          /*
          e[3] = e[2] = e[1] = e[0]
          */
          w2[3][0] = 0.25 * x[3] / (1.0 + x[3] * x[3]);
          w2[3][1] = 0.25 / (1.0 + x[3] * x[3]);
          for (ir = 0; ir < 2; ir++) {
            w2[2][ir] = w2[3][ir];
            w2[1][ir] = w2[3][ir];
            w2[0][ir] = w2[3][ir];
          }
        }
        else {
          /*
          e[3] = e[2] = e[1]
          */
          libtetrabz_polcmplx_1211(x[3], x[0], w2[3]);
          for (ir = 0; ir < 2; ir++) {
            w2[2][ir] = w2[3][ir];
            w2[1][ir] = w2[3][ir];
          }
          libtetrabz_polcmplx_1222(x[0], x[3], w2[0]);
          /*
          # IF(ANY(w2(1:2,1:4) < 0.0)):
          #   WRITE(*,*) ie
          #   WRITE(*,'(100e15.5)') x[0:4]
          #   WRITE(*,'(2e15.5)') w2(1:2,1:4)
          #   STOP "weighting 4=3=2"
          */
        }
      }
      else if (fabs(x[1] - x[0]) < thr) {
        /*
        e[3] = e[2], e[1] = e[0]
        */
        libtetrabz_polcmplx_1221(x[3], x[1], w2[3]);
        for (ir = 0; ir < 2; ir++) w2[2][ir] = w2[3][ir];
        libtetrabz_polcmplx_1221(x[1], x[3], w2[1]);
        for (ir = 0; ir < 2; ir++) w2[0][ir] = w2[1][ir];
        /*
         IF(ANY(w2(1:2,1:4) < 0.0)):
           WRITE(*,*) ie
           WRITE(*,'(100e15.5)') x[0:4]
           WRITE(*,'(2e15.5)') w2(1:2,1:4)
           STOP "weighting 4=3 2=1"
        */
      }
      else {
        /*
        e[3] = e[2]
        */
        libtetrabz_polcmplx_1231(x[3], x[0], x[1], w2[3]);
        for (ir = 0; ir < 2; ir++) w2[2][ir] = w2[3][ir];
        libtetrabz_polcmplx_1233(x[1], x[0], x[3], w2[1]);
        libtetrabz_polcmplx_1233(x[0], x[1], x[3], w2[0]);
        /*
         IF(ANY(w2(1:2,1:4) < 0.0)):
           WRITE(*,*) ie
           WRITE(*,'(100e15.5)') x[0:4]
           WRITE(*,'(2e15.5)') w2(1:2,1:4)
           STOP "weighting 4=3"
        */
      }
    }
    else if (fabs(x[2] - x[1]) < thr) {
      if (fabs(x[2] - x[0]) < thr) {
        /*
        e[2] = e[1] = e[0]
        */
        libtetrabz_polcmplx_1222(x[3], x[2], w2[3]);
        libtetrabz_polcmplx_1211(x[2], x[3], w2[2]);
        for (ir = 0; ir < 2; ir++) w2[1][ir] = w2[2][ir];
        for (ir = 0; ir < 2; ir++) w2[0][ir] = w2[2][ir];
        /*
        IF(ANY(w2(1:2,1:4) < 0.0)):
          WRITE(*,*) ie
          WRITE(*,'(100e15.5)') x[0:4]
          WRITE(*,'(2e15.5)') w2(1:2,1:4)
          STOP "weighting 3=2=1"
        */
      }
      else {
        /*
        e[2] = e[1]
        */
        libtetrabz_polcmplx_1233(x[3], x[0], x[2], w2[3]);
        libtetrabz_polcmplx_1231(x[2], x[0], x[3], w2[2]);
        for (ir = 0; ir < 2; ir++) w2[1][ir] = w2[2][ir];
        libtetrabz_polcmplx_1233(x[0], x[3], x[2], w2[0]);
        /*
        IF(ANY(w2(1:2,1:4) < 0.0)):
          WRITE(*,*) ie
          WRITE(*,'(100e15.5)') x[0:4]
          WRITE(*,'(2e15.5)') w2(1:2,1:4)
          STOP "weighting 3=2"
        */
      }
    }
    else if (fabs(x[1] - x[0]) < thr) {
      /*
      e[1] = e[0]
      */
      libtetrabz_polcmplx_1233(x[3], x[2], x[1], w2[3]);
      libtetrabz_polcmplx_1233(x[2], x[3], x[1], w2[2]);
      libtetrabz_polcmplx_1231(x[1], x[2], x[3], w2[1]);
      for (ir = 0; ir < 2; ir++) w2[0][ir] = w2[1][ir];
      /*
      IF(ANY(w2(1:2,1:4) < 0.0)):
        WRITE(*,*) ie
        WRITE(*,'(100e15.5)') x[0:4]
        WRITE(*,'(2e15.5)') w2(1:2,1:4)
        STOP "weighting 2=1"
      */
    }
    else {
      /*
      Different each other.
      */
      libtetrabz_polcmplx_1234(x[3], x[0], x[1], x[2], w2[3]);
      libtetrabz_polcmplx_1234(x[2], x[0], x[1], x[3], w2[2]);
      libtetrabz_polcmplx_1234(x[1], x[0], x[2], x[3], w2[1]);
      libtetrabz_polcmplx_1234(x[0], x[1], x[2], x[3], w2[0]);
      /*
      IF(ANY(w2(1:2,1:4) < 0.0)):
        WRITE(*,*) ie
        WRITE(*,'(100e15.5)') x[0:4]
        WRITE(*,'(2e15.5)') w2(1:2,1:4)
        STOP "weighting"
      */
    }
    for (i4 = 0; i4 < 4; i4++) {
      w1[indx[i4]][ie][0] = w2[i4][0] / e0[ie][1];
      w1[indx[i4]][ie][1] = w2[i4][1] / (-e0[ie][1]);
    }
  }
}
/*
Tetrahedron method for theta( - E2)
*/
static void libtetrabz_polcmplx2(
  int nb,
  int ne,
  double* *e0,
  double* ei1,
  double** ej1,
  double**** w1
) {
  int ib, i4, j4, ir, ie, indx[4];
  double e[4], tsmall[4][4], v, de[4], thr, *** w2;
  
  w2 = (double***)malloc(4 * sizeof(double**));
  w2[0] = (double**)malloc(4 * ne * sizeof(double*));
  w2[0][0] = (double*)malloc(4 * ne * 2 * sizeof(double));
  for (i4 = 0; i4 < 4; i4++) {
    w2[i4] = w2[0] + i4 * ne;
    for (ie = 0; ie < ne; ie++) {
      w2[i4][ie] = w2[0][0] + i4 * ne * 2 + ie * 2;
    }
  }

  thr = 1.0e-8;

  for (ib = 0; ib < nb; ib++) {

    for (i4 = 0; i4 < 4; i4++)
      for (ie = 0; ie < ne; ie++)
        for (ir = 0; ir < 2; ir++)
          w1[ib][i4][ie][ir] = 0.0;

    for (i4 = 0; i4 < 4; i4++) e[i4] = -ej1[ib][i4];
    eig_sort(4, e, indx);

    if ((e[0] <= 0.0 && 0.0 < e[1]) || (e[0] < 0.0 && 0.0 <= e[1])) {

      libtetrabz_tsmall_a1(e, 0.0, &v, tsmall);

      if (v > thr) {
        for (j4 = 0; j4 < 4; j4++) de[j4] = 0.0;
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            de[j4] += (ej1[ib][indx[i4]] - ei1[indx[i4]]) * tsmall[i4][j4];
        libtetrabz_polcmplx3(ne, e0, de, w2);
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            for (ie = 0; ie < ne; ie++)
              for (ir = 0; ir < 2; ir++)
                w1[ib][indx[i4]][ie][ir] += v * tsmall[i4][j4] * w2[j4][ie][ir];
      }
    }
    else if ((e[1] <= 0.0 && 0.0 < e[2]) || (e[1] < 0.0 && 0.0 <= e[2])) {

      libtetrabz_tsmall_b1(e, 0.0, &v, tsmall);

      if (v > thr) {
        for (j4 = 0; j4 < 4; j4++) de[j4] = 0.0;
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            de[j4] += (ej1[ib][indx[i4]] - ei1[indx[i4]]) * tsmall[i4][j4];
        libtetrabz_polcmplx3(ne, e0, de, w2);
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            for (ie = 0; ie < ne; ie++)
              for (ir = 0; ir < 2; ir++)
                w1[ib][indx[i4]][ie][ir] += v * tsmall[i4][j4] * w2[j4][ie][ir];
      }

      libtetrabz_tsmall_b2(e, 0.0, &v, tsmall);

      if (v > thr) {
        for (j4 = 0; j4 < 4; j4++) de[j4] = 0.0;
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            de[j4] += (ej1[ib][indx[i4]] - ei1[indx[i4]]) * tsmall[i4][j4];
        libtetrabz_polcmplx3(ne, e0, de, w2);
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            for (ie = 0; ie < ne; ie++)
              for (ir = 0; ir < 2; ir++)
                w1[ib][indx[i4]][ie][ir] += v * tsmall[i4][j4] * w2[j4][ie][ir];
      }

      libtetrabz_tsmall_b3(e, 0.0, &v, tsmall);

      if (v > thr) {
        for (j4 = 0; j4 < 4; j4++) de[j4] = 0.0;
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            de[j4] += (ej1[ib][indx[i4]] - ei1[indx[i4]]) * tsmall[i4][j4];
        libtetrabz_polcmplx3(ne, e0, de, w2);
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            for (ie = 0; ie < ne; ie++)
              for (ir = 0; ir < 2; ir++)
                w1[ib][indx[i4]][ie][ir] += v * tsmall[i4][j4] * w2[j4][ie][ir];
      }
    }
    else if ((e[2] <= 0.0 && 0.0 < e[3]) || (e[2] < 0.0 && 0.0 <= e[3])) {

      libtetrabz_tsmall_c1(e, 0.0, &v, tsmall);

      if (v > thr) {
        for (j4 = 0; j4 < 4; j4++) de[j4] = 0.0;
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            de[j4] += (ej1[ib][indx[i4]] - ei1[indx[i4]]) * tsmall[i4][j4];
        libtetrabz_polcmplx3(ne, e0, de, w2);
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            for (ie = 0; ie < ne; ie++)
              for (ir = 0; ir < 2; ir++)
                w1[ib][indx[i4]][ie][ir] += v * tsmall[i4][j4] * w2[j4][ie][ir];
      }

      libtetrabz_tsmall_c2(e, 0.0, &v, tsmall);

      if (v > thr) {
        for (j4 = 0; j4 < 4; j4++) de[j4] = 0.0;
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            de[j4] += (ej1[ib][indx[i4]] - ei1[indx[i4]]) * tsmall[i4][j4];
        libtetrabz_polcmplx3(ne, e0, de, w2);
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            for (ie = 0; ie < ne; ie++)
              for (ir = 0; ir < 2; ir++)
                w1[ib][indx[i4]][ie][ir] += v * tsmall[i4][j4] * w2[j4][ie][ir];
      }

      libtetrabz_tsmall_c3(e, 0.0, &v, tsmall);

      if (v > thr) {
        for (j4 = 0; j4 < 4; j4++) de[j4] = 0.0;
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            de[j4] += (ej1[ib][indx[i4]] - ei1[indx[i4]]) * tsmall[i4][j4];
        libtetrabz_polcmplx3(ne, e0, de, w2);
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            for (ie = 0; ie < ne; ie++)
              for (ir = 0; ir < 2; ir++)
                w1[ib][indx[i4]][ie][ir] += v * tsmall[i4][j4] * w2[j4][ie][ir];
      }
    }
    else if (e[3] <= 0.0) {
      for (i4 = 0; i4 < 4; i4++)
        de[i4] = ej1[ib][i4] - ei1[i4];
      libtetrabz_polcmplx3(ne, e0, de, w2);
      for (i4 = 0; i4 < 4; i4++)
        for (ie = 0; ie < ne; ie++)
          for (ir = 0; ir < 2; ir++)
            w1[ib][i4][ie][ir] += w2[i4][ie][ir];
    }
  }
  free(w2[0][0]);
  free(w2[0]);
  free(w2);
}
/*
Main SUBROUTINE for Polarization (Imaginary axis) : Theta(- E1) * Theta(E2) / (E2 - E1 - iw)
*/
static PyObject* polcmplx_c(PyObject* self, PyObject* args) {
  int it, ik, ie, ib, i4, j4, ir, jb, **ikv, indx[4], i20, ierr, ng[3], nk, nb, ne;
  double wlsm[20][4], **ei1, **ej1, ei2[4], ** ej2, e[4], ***** w1, **** w2, v, tsmall[4][4], thr,
    bvec[3][3], ** eig1, ** eig2, ** e0, ***** wght;
  PyObject* eig1_po, * eig2_po, * e0_po, * wght_po;
  /*
  Read input from python object
  */
  if (!PyArg_ParseTuple(args, "iiiiiidddddddddOOO",
    &ng[0], &ng[1], &ng[2], &nk, &nb, &ne, 
    &bvec[0][0], &bvec[0][1], &bvec[0][2], &bvec[1][0], &bvec[1][1], &bvec[1][2], &bvec[2][0], &bvec[2][1], &bvec[2][2],
    &eig1_po, &eig2_po, &e0_po))
    return NULL;
  /*
  convert python list object to array
  */
  eig1 = (double**)malloc(nk * sizeof(double*));
  eig1[0] = (double*)malloc(nk * nb * sizeof(double));
  eig2 = (double**)malloc(nk * sizeof(double*));
  eig2[0] = (double*)malloc(nk * nb * sizeof(double));
  wght = (double*****)malloc(nk * sizeof(double****));
  wght[0] = (double****)malloc(nk * nb * sizeof(double***));
  wght[0][0] = (double***)malloc(nk * nb * nb * sizeof(double**));
  wght[0][0][0] = (double**)malloc(nk * nb * nb * ne * sizeof(double*));
  wght[0][0][0][0] = (double*)malloc(nk * nb * nb * ne * 2 * sizeof(double));
  for (ik = 0; ik < nk; ik++) {
    eig1[ik] = eig1[0] + ik * nb;
    eig2[ik] = eig2[0] + ik * nb;
    wght[ik] = wght[0] + ik * nb;
    for (ib = 0; ib < nb; ib++) {
      wght[ik][ib] = wght[0][0] + ik * nb * nb + ib * nb;
      for (jb = 0; jb < nb; jb++) {
        wght[ik][ib][jb] = wght[0][0][0]
          + ik * nb * nb * ne + ib * nb * ne + jb * ne;
        for (ie = 0; ie < ne; ie++) {
          wght[ik][ib][jb][ie] = wght[0][0][0][0]
            + ik * nb * nb * ne * 2 + ib * nb * ne * 2 + jb * ne * 2 + ie * 2;
        }
      }
    }
  }
  e0 = (double**)malloc(ne * sizeof(double*));
  e0[0] = (double*)malloc(ne * 2 * sizeof(double));
  for (ie = 0; ie < ne; ie++) 
    e0[ie] = e0[0] + ie * 2;

  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++) {
      eig1[ik][ib] = PyFloat_AsDouble(PyList_GetItem(eig1_po, ik * nb + ib));
      eig2[ik][ib] = PyFloat_AsDouble(PyList_GetItem(eig2_po, ik * nb + ib));
    }
  for (ie = 0; ie < ne; ie++)
    for (ir = 0; ir < 2; ir++)
      e0[ie][ir] = PyFloat_AsDouble(PyList_GetItem(e0_po, ie * 2 + ir));
  /*
  Start main calculation
  */
  ikv = (int**)malloc(6 * nk * sizeof(int*));
  ikv[0] = (int*)malloc(6 * nk * 20 * sizeof(int));
  for (it = 0; it < 6 * nk; it++) {
    ikv[it] = ikv[0] + it * 20;
  }

  ei1 = (double**)malloc(4 * sizeof(double*));
  ej1 = (double**)malloc(4 * sizeof(double*));
  ei1[0] = (double*)malloc(4 * nb * sizeof(double));
  ej1[0] = (double*)malloc(4 * nb * sizeof(double));
  for (i4 = 0; i4 < 4; i4++) {
    ei1[i4] = ei1[0] + i4 * nb;
    ej1[i4] = ej1[0] + i4 * nb;
  }

  ej2 = (double**)malloc(nb * sizeof(double*));
  ej2[0] = (double*)malloc(nb * 4 * sizeof(double));
  for (ib = 0; ib < nb; ib++)
    ej2[ib] = ej2[0] + ib * 4;

  w1 = (double*****)malloc(nb * sizeof(double****));
  w1[0] = (double****)malloc(nb * 4 * sizeof(double***));
  w1[0][0] = (double***)malloc(nb * 4 * nb * sizeof(double**));
  w1[0][0][0] = (double**)malloc(nb * 4 * nb * ne * sizeof(double*));
  w1[0][0][0][0] = (double*)malloc(nb * 4 * nb * ne * 2 * sizeof(double));
  for (ib = 0; ib < nb; ib++) {
    w1[ib] = w1[0] + ib * 4;
    for (i4 = 0; i4 < 4; i4++) {
      w1[ib][i4] = w1[0][0] + ib * 4 * nb + i4 * nb;
      for (jb = 0; jb < nb; jb++) {
        w1[ib][i4][jb] = w1[0][0][0] + ib * 4 * nb * ne + i4 * nb * ne + jb * ne;
        for (ie = 0; ie < ne; ie++) {
          w1[ib][i4][jb][ie] = w1[0][0][0][0] + ib * 4 * nb * ne * 2 + i4 * nb * ne * 2 + jb * ne * 2 + ie * 2;
        }
      }
    }
  }

  w2 = (double****)malloc(nb * sizeof(double***));
  w2[0] = (double***)malloc(nb * 4 * sizeof(double**));
  w2[0][0] = (double**)malloc(nb * 4 * ne * sizeof(double*));
  w2[0][0][0] = (double*)malloc(nb * 4 * ne * 2 * sizeof(double));
  for (ib = 0; ib < nb; ib++) {
    w2[ib] = w2[0] + ib * 4;
    for (i4 = 0; i4 < 4; i4++) {
      w2[ib][i4] = w2[0][0] + ib * 4 * ne + i4 * ne;
      for (ie = 0; ie < ne; ie++) {
        w2[ib][i4][ie] = w2[0][0][0] + ib * 4 * ne * 2 + i4 * ne * 2 + ie * 2;
      }
    }
  }

  libtetrabz_initialize(ng, bvec, wlsm, ikv);

  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      for (jb = 0; jb < nb; jb++)
        for (ie = 0; ie < ne; ie++)
          for (ir = 0; ir < 2; ir++)
            wght[ik][ib][jb][ie][ir] = 0.0;

  thr = 1.0e-8;

  for (it = 0; it < 6 * nk; it++) {

    for (i4 = 0; i4 < 4; i4++)
      for (ib = 0; ib < nb; ib++) {
        ei1[i4][ib] = 0.0;
        ej1[i4][ib] = 0.0;
      }
    for (i20 = 0; i20 < 20; i20++) {
      for (i4 = 0; i4 < 4; i4++) {
        for (ib = 0; ib < nb; ib++) {
          ei1[i4][ib] += eig1[ikv[it][i20]][ib] * wlsm[i20][i4];
          ej1[i4][ib] += eig2[ikv[it][i20]][ib] * wlsm[i20][i4];
        }
      }
    }

    for (ib = 0; ib < nb; ib++) {

      for (i4 = 0; i4 < 4; i4++)
        for (jb = 0; jb < nb; jb++)
          for (ie = 0; ie < ne; ie++)
            for (ir = 0; ir < 2; ir++)
              w1[ib][i4][jb][ie][ir] = 0.0;

      for (i4 = 0; i4 < 4; i4++) e[i4] = ei1[i4][ib];
      eig_sort(4, e, indx);

      if (e[0] <= 0.0 && 0.0 < e[1]) {

        libtetrabz_tsmall_a1(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (j4 = 0; j4 < 4; j4++) {
            ei2[j4] = 0.0;
            for (jb = 0; jb < nb; jb++) ej2[jb][j4] = 0.0;
          }
          for (i4 = 0; i4 < 4; i4++)
            for (j4 = 0; j4 < 4; j4++) {
              ei2[j4] += ei1[indx[i4]][ib] * tsmall[i4][j4];
              for (jb = 0; jb < nb; jb++)
                ej2[jb][j4] += ej1[indx[i4]][jb] * tsmall[i4][j4];
            }
          libtetrabz_polcmplx2(nb, ne, e0, ei2, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j4 = 0; j4 < 4; j4++)
                for (ie = 0; ie < ne; ie++)
                  for (ir = 0; ir < 2; ir++)
                    w1[ib][indx[i4]][jb][ie][ir] += v * tsmall[i4][j4] * w2[jb][j4][ie][ir];
        }
      }
      else if (e[1] <= 0.0 && 0.0 < e[2]) {

        libtetrabz_tsmall_b1(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (j4 = 0; j4 < 4; j4++) {
            ei2[j4] = 0.0;
            for (jb = 0; jb < nb; jb++) ej2[jb][j4] = 0.0;
          }
          for (i4 = 0; i4 < 4; i4++)
            for (j4 = 0; j4 < 4; j4++) {
              ei2[j4] += ei1[indx[i4]][ib] * tsmall[i4][j4];
              for (jb = 0; jb < nb; jb++)
                ej2[jb][j4] += ej1[indx[i4]][jb] * tsmall[i4][j4];
            }
          libtetrabz_polcmplx2(nb, ne, e0, ei2, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j4 = 0; j4 < 4; j4++)
                for (ie = 0; ie < ne; ie++)
                  for (ir = 0; ir < 2; ir++)
                    w1[ib][indx[i4]][jb][ie][ir] += v * tsmall[i4][j4] * w2[jb][j4][ie][ir];
        }

        libtetrabz_tsmall_b2(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (j4 = 0; j4 < 4; j4++) {
            ei2[j4] = 0.0;
            for (jb = 0; jb < nb; jb++) ej2[jb][j4] = 0.0;
          }
          for (i4 = 0; i4 < 4; i4++)
            for (j4 = 0; j4 < 4; j4++) {
              ei2[j4] += ei1[indx[i4]][ib] * tsmall[i4][j4];
              for (jb = 0; jb < nb; jb++)
                ej2[jb][j4] += ej1[indx[i4]][jb] * tsmall[i4][j4];
            }
          libtetrabz_polcmplx2(nb, ne, e0, ei2, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j4 = 0; j4 < 4; j4++)
                for (ie = 0; ie < ne; ie++)
                  for (ir = 0; ir < 2; ir++)
                    w1[ib][indx[i4]][jb][ie][ir] += v * tsmall[i4][j4] * w2[jb][j4][ie][ir];
        }

        libtetrabz_tsmall_b3(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (j4 = 0; j4 < 4; j4++) {
            ei2[j4] = 0.0;
            for (jb = 0; jb < nb; jb++) ej2[jb][j4] = 0.0;
          }
          for (i4 = 0; i4 < 4; i4++)
            for (j4 = 0; j4 < 4; j4++) {
              ei2[j4] += ei1[indx[i4]][ib] * tsmall[i4][j4];
              for (jb = 0; jb < nb; jb++)
                ej2[jb][j4] += ej1[indx[i4]][jb] * tsmall[i4][j4];
            }
          libtetrabz_polcmplx2(nb, ne, e0, ei2, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j4 = 0; j4 < 4; j4++)
                for (ie = 0; ie < ne; ie++)
                  for (ir = 0; ir < 2; ir++)
                    w1[ib][indx[i4]][jb][ie][ir] += v * tsmall[i4][j4] * w2[jb][j4][ie][ir];
        }
      }
      else if (e[2] <= 0.0 && 0.0 < e[3]) {

        libtetrabz_tsmall_c1(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (j4 = 0; j4 < 4; j4++) {
            ei2[j4] = 0.0;
            for (jb = 0; jb < nb; jb++) ej2[jb][j4] = 0.0;
          }
          for (i4 = 0; i4 < 4; i4++)
            for (j4 = 0; j4 < 4; j4++) {
              ei2[j4] += ei1[indx[i4]][ib] * tsmall[i4][j4];
              for (jb = 0; jb < nb; jb++)
                ej2[jb][j4] += ej1[indx[i4]][jb] * tsmall[i4][j4];
            }
          libtetrabz_polcmplx2(nb, ne, e0, ei2, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j4 = 0; j4 < 4; j4++)
                for (ie = 0; ie < ne; ie++)
                  for (ir = 0; ir < 2; ir++)
                    w1[ib][indx[i4]][jb][ie][ir] += v * tsmall[i4][j4] * w2[jb][j4][ie][ir];
        }

        libtetrabz_tsmall_c2(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (j4 = 0; j4 < 4; j4++) {
            ei2[j4] = 0.0;
            for (jb = 0; jb < nb; jb++) ej2[jb][j4] = 0.0;
          }
          for (i4 = 0; i4 < 4; i4++)
            for (j4 = 0; j4 < 4; j4++) {
              ei2[j4] += ei1[indx[i4]][ib] * tsmall[i4][j4];
              for (jb = 0; jb < nb; jb++)
                ej2[jb][j4] += ej1[indx[i4]][jb] * tsmall[i4][j4];
            }
          libtetrabz_polcmplx2(nb, ne, e0, ei2, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j4 = 0; j4 < 4; j4++)
                for (ie = 0; ie < ne; ie++)
                  for (ir = 0; ir < 2; ir++)
                    w1[ib][indx[i4]][jb][ie][ir] += v * tsmall[i4][j4] * w2[jb][j4][ie][ir];
        }

        libtetrabz_tsmall_c3(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (j4 = 0; j4 < 4; j4++) {
            ei2[j4] = 0.0;
            for (jb = 0; jb < nb; jb++) ej2[jb][j4] = 0.0;
          }
          for (i4 = 0; i4 < 4; i4++)
            for (j4 = 0; j4 < 4; j4++) {
              ei2[j4] += ei1[indx[i4]][ib] * tsmall[i4][j4];
              for (jb = 0; jb < nb; jb++)
                ej2[jb][j4] += ej1[indx[i4]][jb] * tsmall[i4][j4];
            }
          libtetrabz_polcmplx2(nb, ne, e0, ei2, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j4 = 0; j4 < 4; j4++)
                for (ie = 0; ie < ne; ie++)
                  for (ir = 0; ir < 2; ir++)
                    w1[ib][indx[i4]][jb][ie][ir] += v * tsmall[i4][j4] * w2[jb][j4][ie][ir];
        }
      }
      else if (e[3] <= 0.0) {
        for (i4 = 0; i4 < 4; i4++) {
          ei2[i4] = ei1[i4][ib];
          for (jb = 0; jb < nb; jb++)
            ej2[jb][i4] = ej1[i4][jb];
        }
        libtetrabz_polcmplx2(nb, ne, e0, ei2, ej2, w2);
        for (i4 = 0; i4 < 4; i4++)
          for (jb = 0; jb < nb; jb++)
            for (ie = 0; ie < ne; ie++)
              for (ir = 0; ir < 2; ir++)
                w1[ib][i4][jb][ie][ir] += w2[jb][i4][ie][ir];
      }
      else {
        continue;
      }
    }
    for (i20 = 0; i20 < 20; i20++)
      for (ib = 0; ib < nb; ib++)
        for (i4 = 0; i4 < 4; i4++)
          for (jb = 0; jb < nb; jb++)
            for (ie = 0; ie < ne; ie++)
              for (ir = 0; ir < 2; ir++)
                wght[ikv[it][i20]][ib][jb][ie][ir] += wlsm[i20][i4] * w1[ib][i4][jb][ie][ir];
  }
  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      for (jb = 0; jb < nb; jb++)
        for (ie = 0; ie < ne; ie++)
          for (i4 = 0; i4 < 2; i4++)
            wght[ik][ib][jb][ie][i4] /= (6.0 * (double) nk);

  free(ikv[0]);
  free(ikv);
  free(ei1[0]);
  free(ei1);
  free(ej1[0]);
  free(ej1);
  free(ej2[0]);
  free(ej2);
  free(w1[0][0][0][0]);
  free(w1[0][0][0]);
  free(w1[0][0]);
  free(w1[0]);
  free(w1);
  free(w2[0][0][0]);
  free(w2[0][0]);
  free(w2[0]);
  free(w2);
  /*
  Convert weight to python list object
  */
  wght_po = PyList_New(nk * nb * nb * ne * 2);
  ierr = 0;

  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      for (jb = 0; jb < nb; jb++)
        for (ie = 0; ie < ne; ie++)
          for (ir = 0; ir < 2; ir++)
            ierr = PyList_SetItem(wght_po, 
              ik * nb * nb * ne * 2 + ib * nb * ne * 2 + jb * ne * 2 + ie * 2 + ir,
              PyFloat_FromDouble(wght[ik][ib][jb][ie][ir]));
  if (ierr != 0)printf("Error in PyList_SetItem\n");

  free(eig1[0]);
  free(eig1);
  free(eig2[0]);
  free(eig2);
  free(wght[0][0][0][0]);
  free(wght[0][0][0]);
  free(wght[0][0]);
  free(wght[0]);
  free(wght);
  free(e0[0]);
  free(e0);

  return wght_po;
}
/*
 Results of Integration (1-x-y-z)/(g0+(g1-g0)x+(g2-g0)y+(g3-g0))
 for 0<x<1, 0<y<1-x, 0<z<1-x-y
*/
/*
1, Different each other
*/
static double libtetrabz_polstat_1234(
    double g1,
    double g2,
    double g3,
    double g4,
    double lng1,
    double lng2,
    double lng3,
    double lng4)
{
  double w2, w3, w4, w;
  w2 = ((lng2 - lng1) / (g2 - g1) * g2 - 1.0) * g2 / (g2 - g1);
  w3 = ((lng3 - lng1) / (g3 - g1) * g3 - 1.0) * g3 / (g3 - g1);
  w4 = ((lng4 - lng1) / (g4 - g1) * g4 - 1.0) * g4 / (g4 - g1);
  w2 = ((w2 - w3) * g2) / (g2 - g3);
  w4 = ((w4 - w3) * g4) / (g4 - g3);
  w = (w4 - w2) / (g4 - g2);
  return w;
}
/*
2, g4 = g1
*/
static double libtetrabz_polstat_1231(
    double g1,
    double g2,
    double g3,
    double lng1,
    double lng2,
    double lng3
) {
  double w2, w3, w;
  w2 = ((lng2 - lng1) / (g2 - g1) * g2 - 1.0) * g2 * g2 / (g2 - g1) - g1 / 2.0;
  w2 = w2 / (g2 - g1);
  w3 = ((lng3 - lng1) / (g3 - g1) * g3 - 1.0) * g3 * g3 / (g3 - g1) - g1 / 2.0;
  w3 = w3 / (g3 - g1);
  w = (w3 - w2) / (g3 - g2);
  return w;
}
/*
# 3, g4 = g3
*/
static double libtetrabz_polstat_1233(
    double g1,
    double g2,
    double g3,
    double lng1,
    double lng2,
    double lng3
) {
  double w2, w3, w;
  w2 = (lng2 - lng1) / (g2 - g1) * g2 - 1.0;
  w2 = (g2 * w2) / (g2 - g1);
  w3 = (lng3 - lng1) / (g3 - g1) * g3 - 1.0;
  w3 = (g3 * w3) / (g3 - g1);
  w2 = (w3 - w2) / (g3 - g2);
  w3 = (lng3 - lng1) / (g3 - g1) * g3 - 1.0;
  w3 = 1.0 - (2.0 * w3 * g1) / (g3 - g1);
  w3 = w3 / (g3 - g1);
  w = (g3 * w3 - g2 * w2) / (g3 - g2);
  return w;
}
/*
4, g4 = g1 and g3 = g2
*/
static double libtetrabz_polstat_1221(
    double g1,
    double g2,
    double lng1,
    double lng2
) {
  double w;
  w = 1.0 - (lng2 - lng1) / (g2 - g1) * g1;
  w = -1.0 + (2.0 * g2 * w) / (g2 - g1);
  w = -1.0 + (3.0 * g2 * w) / (g2 - g1);
  w = w / (2.0 * (g2 - g1));
  return w;
}
/*
5, g4 = g3 = g2
*/
static double libtetrabz_polstat_1222(
    double g1,
    double g2,
    double lng1,
    double lng2
) {
  double w;
  w = (lng2 - lng1) / (g2 - g1) * g2 - 1.0;
  w = (2.0 * g1 * w) / (g2 - g1) - 1.0;
  w = (3.0 * g1 * w) / (g2 - g1) + 1.0;
  w = w / (2.0 * (g2 - g1));
  return w;
}
/*
6, g4 = g3 = g1
*/
static double libtetrabz_polstat_1211(
    double g1,
    double g2,
    double lng1,
    double lng2
) {
  double w;
  w = -1.0 + (lng2 - lng1) / (g2 - g1) * g2;
  w = -1.0 + (2.0 * g2 * w) / (g2 - g1);
  w = -1.0 + (3.0 * g2 * w) / (2.0 * (g2 - g1));
  w = w / (3.0 * (g2 - g1));
  return w;
}
/*
Tetrahedron method for delta(om - ep + e)
*/
static void libtetrabz_polstat3(
    double de[4],
    double w1[4]
)
{
  int i4, indx[4];
  double thr, thr2, e[4], ln[4];

  for (i4 = 0; i4 < 4; i4++) e[i4] = de[i4];
  eig_sort(4, e, indx);

  thr = e[3] * 1.0e-3;
  thr2 = 1.0e-8;

  for(i4 =0; i4 <4; i4++){
    if (e[i4] < thr2) {
      if (i4 == 3) {
        printf("  Nesting # \n");
      }
      ln[i4] = 0.0;
      e[i4] = 0.0;
    }
    else{
      ln[i4] = log(e[i4]);
    }
  }

  if (fabs(e[3] - e[2]) < thr) {
    if (fabs(e[3] - e[1]) < thr) {
      if (fabs(e[3] - e[0]) < thr) {
        /*
        e[3] = e[2] = e[1] = e[0]
        */
        w1[indx[3]] = 0.25 / e[3];
        w1[indx[2]] = w1[indx[3]];
        w1[indx[1]] = w1[indx[3]];
        w1[indx[0]] = w1[indx[3]];
      }
      else {
        /*
        e[3] = e[2] = e[1]
        */
        w1[indx[3]] = libtetrabz_polstat_1211(e[3], e[0], ln[3], ln[0]);
        w1[indx[2]] = w1[indx[3]];
        w1[indx[1]] = w1[indx[3]];
        w1[indx[0]] = libtetrabz_polstat_1222(e[0], e[3], ln[0], ln[3]);
        if (w1[0] < 0.0 || w1[1] < 0.0 || w1[2] < 0.0 || w1[3] < 0.0) {
          printf("%f %f %f %f\n", e[0], e[1], e[2], e[3]);
          printf("%f %f %f %f\n", w1[indx[0]], w1[indx[1]], w1[indx[2]], w1[indx[3]]);
          printf("weighting 4=3=2\n");
        }
      }
    }
    else if (fabs(e[1] - e[0]) < thr) {
      /*
      e[3] = e[2], e[1] = e[0]
      */
      w1[indx[3]] = libtetrabz_polstat_1221(e[3], e[1], ln[3], ln[1]);
      w1[indx[2]] = w1[indx[3]];
      w1[indx[1]] = libtetrabz_polstat_1221(e[1], e[3], ln[1], ln[3]);
      w1[indx[0]] = w1[indx[1]];

      if (w1[0] < 0.0 || w1[1] < 0.0 || w1[2] < 0.0 || w1[3] < 0.0) {
        printf("%f %f %f %f\n", e[0], e[1], e[2], e[3]);
        printf("%f %f %f %f\n", w1[indx[0]], w1[indx[1]], w1[indx[2]], w1[indx[3]]);
        printf("weighting 4=3 2=1\n");
      }
    }
    else {
      /*
      e[3] = e[2]
      */
      w1[indx[3]] = libtetrabz_polstat_1231(e[3], e[0], e[1], ln[3], ln[0], ln[1]);
      w1[indx[2]] = w1[indx[3]];
      w1[indx[1]] = libtetrabz_polstat_1233(e[1], e[0], e[3], ln[1], ln[0], ln[3]);
      w1[indx[0]] = libtetrabz_polstat_1233(e[0], e[1], e[3], ln[0], ln[1], ln[3]);

      if (w1[0] < 0.0 || w1[1] < 0.0 || w1[2] < 0.0 || w1[3] < 0.0) {
        printf("%f %f %f %f\n", e[0], e[1], e[2], e[3]);
        printf("%f %f %f %f\n", w1[indx[0]], w1[indx[1]], w1[indx[2]], w1[indx[3]]);
        printf("weighting 4=3\n");
      }
    }
  }  
  else if (fabs(e[2] - e[1]) < thr) {
    if (fabs(e[2] - e[0]) < thr) {
      /*
      e[2] = e[1] = e[0]
      */
      w1[indx[3]] = libtetrabz_polstat_1222(e[3], e[2], ln[3], ln[2]);
      w1[indx[2]] = libtetrabz_polstat_1211(e[2], e[3], ln[2], ln[3]);
      w1[indx[1]] = w1[indx[2]];
      w1[indx[0]] = w1[indx[2]];

      if (w1[0] < 0.0 || w1[1] < 0.0 || w1[2] < 0.0 || w1[3] < 0.0) {
        printf("%f %f %f %f\n", e[0], e[1], e[2], e[3]);
        printf("%f %f %f %f\n", w1[indx[0]], w1[indx[1]], w1[indx[2]], w1[indx[3]]);
        printf("weighting 3=2=1\n");
      }
    }
    else {
      /*
      e[2] = e[1]
      */
      w1[indx[3]] = libtetrabz_polstat_1233(e[3], e[0], e[2], ln[3], ln[0], ln[2]);
      w1[indx[2]] = libtetrabz_polstat_1231(e[2], e[0], e[3], ln[2], ln[0], ln[3]);
      w1[indx[1]] = w1[indx[2]];
      w1[indx[0]] = libtetrabz_polstat_1233(e[0], e[3], e[2], ln[0], ln[3], ln[2]);

      if (w1[0] < 0.0 || w1[1] < 0.0 || w1[2] < 0.0 || w1[3] < 0.0) {
        printf("%f %f %f %f\n", e[0], e[1], e[2], e[3]);
        printf("%f %f %f %f\n", w1[indx[0]], w1[indx[1]], w1[indx[2]], w1[indx[3]]);
        printf("weighting 3=2\n");
      }
    }
  }  
  else if (fabs(e[1] - e[0]) < thr) {
    /*
    e[1] = e[0]
    */
    w1[indx[3]] = libtetrabz_polstat_1233(e[3], e[2], e[1], ln[3], ln[2], ln[1]);
    w1[indx[2]] = libtetrabz_polstat_1233(e[2], e[3], e[1], ln[2], ln[3], ln[1]);
    w1[indx[1]] = libtetrabz_polstat_1231(e[1], e[2], e[3], ln[1], ln[2], ln[3]);
    w1[indx[0]] = w1[indx[1]];

    if (w1[0] < 0.0 || w1[1] < 0.0 || w1[2] < 0.0 || w1[3] < 0.0) {
      printf("%f %f %f %f\n", e[0], e[1], e[2], e[3]);
      printf("%f %f %f %f\n", w1[indx[0]], w1[indx[1]], w1[indx[2]], w1[indx[3]]);
      printf("weighting 2=1\n");
    }
  }  
  else {
    /*
    Different each other.
    */
    w1[indx[3]] = libtetrabz_polstat_1234(e[3], e[0], e[1], e[2], ln[3], ln[0], ln[1], ln[2]);
    w1[indx[2]] = libtetrabz_polstat_1234(e[2], e[0], e[1], e[3], ln[2], ln[0], ln[1], ln[3]);
    w1[indx[1]] = libtetrabz_polstat_1234(e[1], e[0], e[2], e[3], ln[1], ln[0], ln[2], ln[3]);
    w1[indx[0]] = libtetrabz_polstat_1234(e[0], e[1], e[2], e[3], ln[0], ln[1], ln[2], ln[3]);

    if (w1[0] < 0.0 || w1[1] < 0.0 || w1[2] < 0.0 || w1[3] < 0.0) {
      printf("%f %f %f %f\n", e[0], e[1], e[2], e[3]);
      printf("%f %f %f %f\n", w1[indx[0]], w1[indx[1]], w1[indx[2]], w1[indx[3]]);
      printf("weighting\n");
    }
  }
}
/*
Tetrahedron method for theta( - E2)
*/
static void libtetrabz_polstat2(
    int nb,
    double *ei1,
    double **ej1,
    double **w1
) {
  int i4, j4, ib, indx[4];
  double v, thr, e[4], de[4], w2[4], tsmall[4][4];

  thr = 1.0e-8;

  for (ib = 0; ib < nb; ib++) {

    for (i4 = 0; i4 < 4; i4++)
      w1[ib][i4] = 0.0;

    for (i4 = 0; i4 < 4; i4++) e[i4] = -ej1[ib][i4];
    eig_sort(4, e, indx);

    if ((e[0] <= 0.0 && 0.0 < e[1]) || (e[0] < 0.0 && 0.0 <= e[1])) {

      libtetrabz_tsmall_a1(e, 0.0, &v, tsmall);

      if (v > thr) {
        for (j4 = 0; j4 < 4; j4++) de[j4] = 0.0;
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            de[j4] += (ej1[ib][indx[i4]] - ei1[indx[i4]]) * tsmall[i4][j4];
        libtetrabz_polstat3(de, w2);
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            w1[ib][indx[i4]] += v * tsmall[i4][j4] * w2[j4];
      }
    }
    else if ((e[1] <= 0.0 && 0.0 < e[2]) || (e[1] < 0.0 && 0.0 <= e[2])) {

      libtetrabz_tsmall_b1(e, 0.0, &v, tsmall);

      if (v > thr) {
        for (j4 = 0; j4 < 4; j4++) de[j4] = 0.0;
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            de[j4] += (ej1[ib][indx[i4]] - ei1[indx[i4]]) * tsmall[i4][j4];
        libtetrabz_polstat3(de, w2);
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            w1[ib][indx[i4]] += v * tsmall[i4][j4] * w2[j4];
      }

      libtetrabz_tsmall_b2(e, 0.0, &v, tsmall);

      if (v > thr) {
        for (j4 = 0; j4 < 4; j4++) de[j4] = 0.0;
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            de[j4] += (ej1[ib][indx[i4]] - ei1[indx[i4]]) * tsmall[i4][j4];
        libtetrabz_polstat3(de, w2);
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            w1[ib][indx[i4]] += v * tsmall[i4][j4] * w2[j4];
      }

      libtetrabz_tsmall_b3(e, 0.0, &v, tsmall);

      if (v > thr) {
        for (j4 = 0; j4 < 4; j4++) de[j4] = 0.0;
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            de[j4] += (ej1[ib][indx[i4]] - ei1[indx[i4]]) * tsmall[i4][j4];
        libtetrabz_polstat3(de, w2);
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            w1[ib][indx[i4]] += v * tsmall[i4][j4] * w2[j4];
      }
    }
    else if ((e[2] <= 0.0 && 0.0 < e[3]) || (e[2] < 0.0 && 0.0 <= e[3])) {

      libtetrabz_tsmall_c1(e, 0.0, &v, tsmall);
      
      if (v > thr) {
        for (j4 = 0; j4 < 4; j4++) de[j4] = 0.0;
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            de[j4] += (ej1[ib][indx[i4]] - ei1[indx[i4]]) * tsmall[i4][j4];
        libtetrabz_polstat3(de, w2);
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            w1[ib][indx[i4]] += v * tsmall[i4][j4] * w2[j4];
      }

      libtetrabz_tsmall_c2(e, 0.0, &v, tsmall);

      if (v > thr) {
        for (j4 = 0; j4 < 4; j4++) de[j4] = 0.0;
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            de[j4] += (ej1[ib][indx[i4]] - ei1[indx[i4]]) * tsmall[i4][j4];
        libtetrabz_polstat3(de, w2);
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            w1[ib][indx[i4]] += v * tsmall[i4][j4] * w2[j4];
      }

      libtetrabz_tsmall_c3(e, 0.0, &v, tsmall);

      if (v > thr) {
        for (j4 = 0; j4 < 4; j4++) de[j4] = 0.0;
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            de[j4] += (ej1[ib][indx[i4]] - ei1[indx[i4]]) * tsmall[i4][j4];
        libtetrabz_polstat3(de, w2);
        for (i4 = 0; i4 < 4; i4++)
          for (j4 = 0; j4 < 4; j4++)
            w1[ib][indx[i4]] += v * tsmall[i4][j4] * w2[j4];
      }
    }
    else if (e[3] <= 0.0) {
      for (i4 = 0; i4 < 4; i4++)
        de[i4] = ej1[ib][i4] - ei1[i4];
      libtetrabz_polstat3(de, w2);
      for (i4 = 0; i4 < 4; i4++) 
        w1[ib][i4] += w2[i4];
    }
  }
}
/*
Compute Static polarization function
*/
static PyObject* polstat_c(PyObject* self, PyObject* args)
{
  int it, ik, ib, jb, i20, i4, j4, **ikv, indx[4], ierr, ng[3], nk, nb;
  double wlsm[20][4], **ei1, **ej1, ei2[4], ** ej2, e[4], *** w1, ** w2, v, tsmall[4][4], thr,
    bvec[3][3], ** eig1, ** eig2, *** wght;
  PyObject* eig1_po, * eig2_po, * wght_po;
  /*
   Read input from python object
   */
  if (!PyArg_ParseTuple(args, "iiiiidddddddddOO",
    &ng[0], &ng[1], &ng[2], &nk, &nb,
    &bvec[0][0], &bvec[0][1], &bvec[0][2], &bvec[1][0], &bvec[1][1], &bvec[1][2], &bvec[2][0], &bvec[2][1], &bvec[2][2],
    &eig1_po, &eig2_po))
    return NULL;
  /*
  convert python list object to array
  */
  eig1 = (double**)malloc(nk * sizeof(double*));
  eig1[0] = (double*)malloc(nk * nb * sizeof(double));
  eig2 = (double**)malloc(nk * sizeof(double*));
  eig2[0] = (double*)malloc(nk * nb * sizeof(double));
  wght = (double***)malloc(nk * sizeof(double**));
  wght[0] = (double**)malloc(nk * nb * sizeof(double*));
  wght[0][0] = (double*)malloc(nk * nb * nb * sizeof(double));
  for (ik = 0; ik < nk; ik++) {
    eig1[ik] = eig1[0] + ik * nb;
    eig2[ik] = eig2[0] + ik * nb;
    wght[ik] = wght[0] + ik * nb;
    for (ib = 0; ib < nb; ib++) {
      wght[ik][ib] = wght[0][0] + ik * nb * nb + ib * nb;
    }
  }

  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++) {
      eig1[ik][ib] = PyFloat_AsDouble(PyList_GetItem(eig1_po, ik * nb + ib));
      eig2[ik][ib] = PyFloat_AsDouble(PyList_GetItem(eig2_po, ik * nb + ib));
    }
  /*
  Start main calculation
  */


  thr = 1.0e-10;

  ikv = (int**)malloc(6 * nk * sizeof(int*));
  ikv[0] = (int*)malloc(6 * nk * 20 * sizeof(int));
  for (it = 0; it < 6 * nk; it++) {
    ikv[it] = ikv[0] + it * 20;
  }

  ei1 = (double**)malloc(4 * sizeof(double*));
  ej1 = (double**)malloc(4 * sizeof(double*));
  ei1[0] = (double*)malloc(4 * nb * sizeof(double));
  ej1[0] = (double*)malloc(4 * nb * sizeof(double));
  for (i4 = 0; i4 < 4; i4++) {
    ei1[i4] = ei1[0] + i4 * nb;
    ej1[i4] = ej1[0] + i4 * nb;
  }

  ej2 = (double**)malloc(nb * sizeof(double*));
  ej2[0] = (double*)malloc(nb * 4 * sizeof(double));
  for (ib = 0; ib < nb; ib++)
    ej2[ib] = ej2[0] + ib * 4;

  w1 = (double***)malloc(nb * sizeof(double**));
  w1[0] = (double**)malloc(nb * 4 * sizeof(double*));
  w1[0][0] = (double*)malloc(nb * 4 * nb * sizeof(double));
  for (ib = 0; ib < nb; ib++) {
    w1[ib] = w1[0] + ib * 4;
    for (i4 = 0; i4 < 4; i4++) {
      w1[ib][i4] = w1[0][0] + ib * 4 * nb + i4 * nb;
    }
  }

  w2 = (double**)malloc(nb * sizeof(double*));
  w2[0] = (double*)malloc(nb * 4 * sizeof(double));
  for (ib = 0; ib < nb; ib++) {
    w2[ib] = w2[0] + ib * 4;
  }

  libtetrabz_initialize(ng, bvec, wlsm, ikv);

  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      for (jb = 0; jb < nb; jb++)
        wght[ik][ib][jb] = 0.0;

  for(it = 0; it < 6*nk; it++) {

    for (i4 = 0; i4 < 4; i4++)
      for (ib = 0; ib < nb; ib++) {
        ei1[i4][ib] = 0.0;
        ej1[i4][ib] = 0.0;
      }
    for (i20 = 0; i20 < 20; i20++) {
      for (i4 = 0; i4 < 4; i4++) {
        for (ib = 0; ib < nb; ib++) {
          ei1[i4][ib] += eig1[ikv[it][i20]][ib] * wlsm[i20][i4];
          ej1[i4][ib] += eig2[ikv[it][i20]][ib] * wlsm[i20][i4];
        }
      }
    }
    
    for (ib = 0; ib < nb; ib++) {

      for (i4 = 0; i4 < 4; i4++)
        for (jb = 0; jb < nb; jb++)
          w1[ib][i4][jb] = 0.0;

      for (i4 = 0; i4 < 4; i4++) e[i4] = ei1[i4][ib];
      eig_sort(4, e, indx);

      if (e[0] <= 0.0 && 0.0 < e[1]) {

        libtetrabz_tsmall_a1(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (j4 = 0; j4 < 4; j4++) {
            ei2[j4] = 0.0;
            for (jb = 0; jb < nb; jb++) ej2[jb][j4] = 0.0;
          }
          for (i4 = 0; i4 < 4; i4++)
            for (j4 = 0; j4 < 4; j4++) {
              ei2[j4] += ei1[indx[i4]][ib] * tsmall[i4][j4];
              for (jb = 0; jb < nb; jb++)
                ej2[jb][j4] += ej1[indx[i4]][jb] * tsmall[i4][j4];
            }
          libtetrabz_polstat2(nb, ei2, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j4 = 0; j4 < 4; j4++)
                w1[ib][indx[i4]][jb] += v * tsmall[i4][j4] * w2[jb][j4];
        }
      }
      else if (e[1] <= 0.0 && 0.0 < e[2]) {

        libtetrabz_tsmall_b1(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (j4 = 0; j4 < 4; j4++) {
            ei2[j4] = 0.0;
            for (jb = 0; jb < nb; jb++) ej2[jb][j4] = 0.0;
          }
          for (i4 = 0; i4 < 4; i4++)
            for (j4 = 0; j4 < 4; j4++) {
              ei2[j4] += ei1[indx[i4]][ib] * tsmall[i4][j4];
              for (jb = 0; jb < nb; jb++)
                ej2[jb][j4] += ej1[indx[i4]][jb] * tsmall[i4][j4];
            }
          libtetrabz_polstat2(nb, ei2, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j4 = 0; j4 < 4; j4++)
                w1[ib][indx[i4]][jb] += v * tsmall[i4][j4] * w2[jb][j4];
        }

        libtetrabz_tsmall_b2(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (j4 = 0; j4 < 4; j4++) {
            ei2[j4] = 0.0;
            for (jb = 0; jb < nb; jb++) ej2[jb][j4] = 0.0;
          }
          for (i4 = 0; i4 < 4; i4++)
            for (j4 = 0; j4 < 4; j4++) {
              ei2[j4] += ei1[indx[i4]][ib] * tsmall[i4][j4];
              for (jb = 0; jb < nb; jb++)
                ej2[jb][j4] += ej1[indx[i4]][jb] * tsmall[i4][j4];
            }
          libtetrabz_polstat2(nb, ei2, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j4 = 0; j4 < 4; j4++)
                w1[ib][indx[i4]][jb] += v * tsmall[i4][j4] * w2[jb][j4];
        }

        libtetrabz_tsmall_b3(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (j4 = 0; j4 < 4; j4++) {
            ei2[j4] = 0.0;
            for (jb = 0; jb < nb; jb++) ej2[jb][j4] = 0.0;
          }
          for (i4 = 0; i4 < 4; i4++)
            for (j4 = 0; j4 < 4; j4++) {
              ei2[j4] += ei1[indx[i4]][ib] * tsmall[i4][j4];
              for (jb = 0; jb < nb; jb++)
                ej2[jb][j4] += ej1[indx[i4]][jb] * tsmall[i4][j4];
            }
          libtetrabz_polstat2(nb, ei2, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j4 = 0; j4 < 4; j4++)
                w1[ib][indx[i4]][jb] += v * tsmall[i4][j4] * w2[jb][j4];
        }
      }
      else if (e[2] <= 0.0 && 0.0 < e[3]) {

        libtetrabz_tsmall_c1(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (j4 = 0; j4 < 4; j4++) {
            ei2[j4] = 0.0;
            for (jb = 0; jb < nb; jb++) ej2[jb][j4] = 0.0;
          }
          for (i4 = 0; i4 < 4; i4++)
            for (j4 = 0; j4 < 4; j4++) {
              ei2[j4] += ei1[indx[i4]][ib] * tsmall[i4][j4];
              for (jb = 0; jb < nb; jb++)
                ej2[jb][j4] += ej1[indx[i4]][jb] * tsmall[i4][j4];
            }
          libtetrabz_polstat2(nb, ei2, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j4 = 0; j4 < 4; j4++)
                w1[ib][indx[i4]][jb] += v * tsmall[i4][j4] * w2[jb][j4];
        }

        libtetrabz_tsmall_c2(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (j4 = 0; j4 < 4; j4++) {
            ei2[j4] = 0.0;
            for (jb = 0; jb < nb; jb++) ej2[jb][j4] = 0.0;
          }
          for (i4 = 0; i4 < 4; i4++)
            for (j4 = 0; j4 < 4; j4++) {
              ei2[j4] += ei1[indx[i4]][ib] * tsmall[i4][j4];
              for (jb = 0; jb < nb; jb++)
                ej2[jb][j4] += ej1[indx[i4]][jb] * tsmall[i4][j4];
            }
          libtetrabz_polstat2(nb, ei2, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j4 = 0; j4 < 4; j4++)
                w1[ib][indx[i4]][jb] += v * tsmall[i4][j4] * w2[jb][j4];
        }

        libtetrabz_tsmall_c3(e, 0.0, &v, tsmall);

        if (v > thr) {
          for (j4 = 0; j4 < 4; j4++) {
            ei2[j4] = 0.0;
            for (jb = 0; jb < nb; jb++) ej2[jb][j4] = 0.0;
          }
          for (i4 = 0; i4 < 4; i4++)
            for (j4 = 0; j4 < 4; j4++) {
              ei2[j4] += ei1[indx[i4]][ib] * tsmall[i4][j4];
              for (jb = 0; jb < nb; jb++)
                ej2[jb][j4] += ej1[indx[i4]][jb] * tsmall[i4][j4];
            }
          libtetrabz_polstat2(nb, ei2, ej2, w2);
          for (i4 = 0; i4 < 4; i4++)
            for (jb = 0; jb < nb; jb++)
              for (j4 = 0; j4 < 4; j4++)
                w1[ib][indx[i4]][jb] += v * tsmall[i4][j4] * w2[jb][j4];
        }
      }
      else if (e[3] <= 0.0) {
        for (i4 = 0; i4 < 4; i4++) {
          ei2[i4] = ei1[i4][ib];
          for (jb = 0; jb < nb; jb++)
            ej2[jb][i4] = ej1[i4][jb];
        }
        libtetrabz_polstat2(nb, ei2, ej2, w2);
        for (i4 = 0; i4 < 4; i4++)
          for (jb = 0; jb < nb; jb++)
            w1[ib][i4][jb] += w2[jb][i4];
      }
      else {
        continue;
      }
    }
    for (i20 = 0; i20 < 20; i20++)
      for (ib = 0; ib < nb; ib++)
        for (i4 = 0; i4 < 4; i4++)
          for (jb = 0; jb < nb; jb++)
            wght[ikv[it][i20]][ib][jb] += wlsm[i20][i4] * w1[ib][i4][jb];
  }
  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      for (jb = 0; jb < nb; jb++)
        wght[ik][ib][jb] /= (6.0 * (double) nk);

  free(ikv[0]);
  free(ikv);
  free(ei1[0]);
  free(ei1);
  free(ej1[0]);
  free(ej1);
  free(ej2[0]);
  free(ej2);
  free(w1[0][0]);
  free(w1[0]);
  free(w1);
  free(w2[0]);
  free(w2);
  /*
  Convert weight to python list object
  */
  wght_po = PyList_New(nk * nb * nb);
  ierr = 0;
  for (ik = 0; ik < nk; ik++)
    for (ib = 0; ib < nb; ib++)
      for (jb = 0; jb < nb; jb++)
        ierr = PyList_SetItem(wght_po, ik * nb * nb + ib * nb + jb, PyFloat_FromDouble(wght[ik][ib][jb]));
  if (ierr != 0)printf("Error in PyList_SetItem\n");

  free(eig1[0]);
  free(eig1);
  free(eig2[0]);
  free(eig2);
  free(wght[0][0]);
  free(wght[0]);
  free(wght);

  return wght_po;
}
static PyMethodDef LibTetraBZCMethods[] = {
    {"dos_c",  dos_c, METH_VARARGS, "Density of States"},
    {"intdos_c",  intdos_c, METH_VARARGS, "Integrated density of States"},
    {"occ_c",  occ_c, METH_VARARGS, "Occupations"},
    {"fermieng_c",  fermieng_c, METH_VARARGS, "Fermi energy"},
    {"dbldelta_c",  dbldelta_c, METH_VARARGS, "Double delta"},
    {"dblstep_c",  dblstep_c, METH_VARARGS, "Double step function"},
    {"polstat_c",  polstat_c, METH_VARARGS, "Static polarization function"},
    {"fermigr_c",  fermigr_c, METH_VARARGS, "Fermi's golden rule"},
    {"polcmplx_c",  polcmplx_c, METH_VARARGS, "Dynamical polarization function on complex frequency"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef libtetrabzcmodule = {
    PyModuleDef_HEAD_INIT,
    "libtetrabzc",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    LibTetraBZCMethods
};

PyMODINIT_FUNC
PyInit_libtetrabzc(void)
{
  return PyModule_Create(&libtetrabzcmodule);
}
