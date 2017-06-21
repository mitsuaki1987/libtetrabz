#include <complex.h>
#pragma once

void libtetrabz_occ(     int *ltetra0, double *bvec, int *nb0, int *nge,
                         double *eig, int *ngw, double *wght0, int *comm);
void libtetrabz_fermieng(int *ltetra0, double *bvec, int *nb0, int *nge,
                         double *eig, int *ngw, double *wght0, double *ef,
                         double *nelec, int *comm);
void libtetrabz_dos(     int *ltetra0, double *bvec, int *nb0, int *nge,
                         double *eig, int *ngw, double *wght0, int *ne0,
                         double *e0, int *comm);
void libtetrabz_dbldelta(int *ltetra0, double *bvec, int *nb0, int *nge,
                         double *eig1, double *eig2, int *ngw, double *wght0,
                         int *comm);
void libtetrabz_dblstep( int *ltetra0, double *bvec, int *nb0, int *nge,
                         double *eig1, double *eig2, int *ngw, double *wght0,
                         int *comm);
void libtetrabz_polstat( int *ltetra0, double *bvec, int *nb0, int *nge,
                         double *eig1, double *eig2, int *ngw, double *wght0,
                         int *comm);
void libtetrabz_fermigr( int *ltetra0, double *bvec, int *nb0, int *nge,
                         double *eig1, double *eig2, int *ngw, double *wght0,
                         int *ne0, double *e0, int *comm);
void libtetrabz_polcmplx(int *ltetra0, double *bvec, int *nb0, int *nge,
                         double *eig1, double *eig2, int *ngw, double *wght0,
                         int *ne0, double complex *e0, int *comm);
