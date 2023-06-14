#include <complex.h>
#pragma once

void libtetrabz_occ(     const int *ltetra0, const double *bvec, const int *nb0, const int *nge,
                         const double *eig, const int *ngw, double *wght0, const int *comm);
void libtetrabz_fermieng(const int *ltetra0, const double *bvec, const int *nb0, const int *nge,
                         const double *eig, const int *ngw, double *wght0, double *ef,
                         const double *nelec, const int *comm);
void libtetrabz_dos(     const int *ltetra0, const double *bvec, const int *nb0, const int *nge,
                         const double *eig, const int *ngw, double *wght0, const int *ne0,
                         const double *e0, const int *comm);
void libtetrabz_dbldelta(const int *ltetra0, const double *bvec, const int *nb0, const int *nge,
                         const double *eig1, const double *eig2, const int *ngw, double *wght0,
                         const int *comm);
void libtetrabz_dblstep( const int *ltetra0, const double *bvec, const int *nb0, const int *nge,
                         const double *eig1, const double *eig2, const int *ngw, double *wght0,
                         const int *comm);
void libtetrabz_polstat( const int *ltetra0, const double *bvec, const int *nb0, const int *nge,
                         const double *eig1, const double *eig2, const int *ngw, double *wght0,
                         const int *comm);
void libtetrabz_fermigr( const int *ltetra0, const double *bvec, const int *nb0, const int *nge,
                         const double *eig1, const double *eig2, const int *ngw, double *wght0,
                         const int *ne0, const double *e0, const int *comm);
void libtetrabz_polcmplx(const int *ltetra0, const double *bvec, const int *nb0, const int *nge,
                         const double *eig1, const double *eig2, const int *ngw, double complex *wght0,
                         const int *ne0, const double complex *e0, const int *comm);
