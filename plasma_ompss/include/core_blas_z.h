/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @precisions normal z -> s d c
 *
 **/
#ifndef ICL_CORE_BLAS_Z_H
#define ICL_CORE_BLAS_Z_H

#include "plasma_async.h"
#include "plasma_barrier.h"
#include "plasma_descriptor.h"
#include "plasma_types.h"
#include "plasma_workspace.h"
#include "plasma_descriptor.h"

#ifdef __cplusplus
extern "C" {
#endif

#define COMPLEX

/******************************************************************************/
#ifdef COMPLEX
double core_dcabs1(plasma_complex64_t alpha);
#endif

void core_zgemm(plasma_enum_t transa, plasma_enum_t transb,
                int m, int n, int k,
                plasma_complex64_t alpha, const plasma_complex64_t *A, int lda,
                                          const plasma_complex64_t *B, int ldb,
                plasma_complex64_t beta,        plasma_complex64_t *C, int ldc);

int core_zgeqrt(int m, int n, int ib,
                plasma_complex64_t *A, int lda,
                plasma_complex64_t *T, int ldt,
                plasma_complex64_t *tau,
                plasma_complex64_t *work);

void core_zgetrf(plasma_desc_t A, int *ipiv, int ib, int rank, int size,
                 volatile int *max_idx, volatile plasma_complex64_t *max_val,
                 volatile int *info, plasma_barrier_t *barrier);

void core_zherk(plasma_enum_t uplo, plasma_enum_t trans,
                int n, int k,
                double alpha, const plasma_complex64_t *A, int lda,
                double beta,        plasma_complex64_t *C, int ldc);

void core_zlacpy(plasma_enum_t uplo, plasma_enum_t transa,
                 int m, int n,
                 const plasma_complex64_t *A, int lda,
                       plasma_complex64_t *B, int ldb);

void core_zlaset(plasma_enum_t uplo,
                 int m, int n,
                 plasma_complex64_t alpha, plasma_complex64_t beta,
                 plasma_complex64_t *A, int lda);

void core_zgeswp(plasma_enum_t colrow,
                 plasma_desc_t A, int k1, int k2, const int *ipiv, int incx);

int core_zpamm(int op, plasma_enum_t side, plasma_enum_t storev,
               int m, int n, int k, int l,
               const plasma_complex64_t *A1, int lda1,
                     plasma_complex64_t *A2, int lda2,
               const plasma_complex64_t *V,  int ldv,
                     plasma_complex64_t *W,  int ldw);

int core_zparfb(plasma_enum_t side, plasma_enum_t trans, plasma_enum_t direct,
                plasma_enum_t storev,
                int m1, int n1, int m2, int n2, int k, int l,
                      plasma_complex64_t *A1,   int lda1,
                      plasma_complex64_t *A2,   int lda2,
                const plasma_complex64_t *V,    int ldv,
                const plasma_complex64_t *T,    int ldt,
                      plasma_complex64_t *work, int ldwork);

int core_zpemv(plasma_enum_t trans, int storev,
               int m, int n, int l,
               plasma_complex64_t alpha,
               const plasma_complex64_t *A, int lda,
               const plasma_complex64_t *X, int incx,
               plasma_complex64_t beta,
               plasma_complex64_t *Y, int incy,
               plasma_complex64_t *work);

int core_zpotrf(plasma_enum_t uplo,
                int n,
                plasma_complex64_t *A, int lda);

void core_zsyrk(plasma_enum_t uplo, plasma_enum_t trans,
                int n, int k,
                plasma_complex64_t alpha, const plasma_complex64_t *A, int lda,
                plasma_complex64_t beta,        plasma_complex64_t *C, int ldc);

void core_ztrsm(plasma_enum_t side, plasma_enum_t uplo,
                plasma_enum_t transa, plasma_enum_t diag,
                int m, int n,
                plasma_complex64_t alpha, const plasma_complex64_t *A, int lda,
                                                plasma_complex64_t *B, int ldb);

int core_ztsmqr(plasma_enum_t side, plasma_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      plasma_complex64_t *A1,   int lda1,
                      plasma_complex64_t *A2,   int lda2,
                const plasma_complex64_t *V,    int ldv,
                const plasma_complex64_t *T,    int ldt,
                      plasma_complex64_t *work, int ldwork);

int core_ztsqrt(int m, int n, int ib,
                plasma_complex64_t *A1, int lda1,
                plasma_complex64_t *A2, int lda2,
                plasma_complex64_t *T,  int ldt,
                plasma_complex64_t *tau,
                plasma_complex64_t *work);

int core_zunmqr(plasma_enum_t side, plasma_enum_t trans,
                int m, int n, int k, int ib,
                const plasma_complex64_t *A,    int lda,
                const plasma_complex64_t *T,    int ldt,
                      plasma_complex64_t *C,    int ldc,
                      plasma_complex64_t *work, int ldwork);

/******************************************************************************/
void core_ompss_zgemm(
    plasma_enum_t transa, plasma_enum_t transb,
    int m, int n, int k,
    plasma_complex64_t alpha, const plasma_complex64_t *A, int lda,
                              const plasma_complex64_t *B, int ldb,
    plasma_complex64_t beta,        plasma_complex64_t *C, int ldc,
    plasma_sequence_t *sequence, plasma_request_t *request);

void core_ompss_zgeqrt(int m, int n, int ib,
                       plasma_complex64_t *A, int lda,
                       plasma_complex64_t *T, int ldt,
                       plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void core_ompss_zherk(plasma_enum_t uplo, plasma_enum_t trans,
                    int n, int k,
                    double alpha, const plasma_complex64_t *A, int lda,
                    double beta,        plasma_complex64_t *C, int ldc,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void core_ompss_zlacpy(plasma_enum_t uplo, plasma_enum_t transa,
                       int m, int n,
                       const plasma_complex64_t *A, int lda,
                             plasma_complex64_t *B, int ldb,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void core_ompss_zlaset(plasma_enum_t uplo,
                       int mb, int nb,
                       int i, int j,
                       int m, int n,
                       plasma_complex64_t alpha, plasma_complex64_t beta,
                       plasma_complex64_t *A);

void core_ompss_zpotrf(plasma_enum_t uplo,
                     int n,
                     plasma_complex64_t *A, int lda,
                     int iinfo,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_ompss_ztrsm(
    plasma_enum_t side, plasma_enum_t uplo,
    plasma_enum_t transa, plasma_enum_t diag,
    int m, int n,
    plasma_complex64_t alpha, const plasma_complex64_t *A, int lda,
                                    plasma_complex64_t *B, int ldb,
    plasma_sequence_t *sequence, plasma_request_t *request);

void core_ompss_ztsmqr(plasma_enum_t side, plasma_enum_t trans,
                       int m1, int n1, int m2, int n2, int k, int ib,
                             plasma_complex64_t *A1, int lda1,
                             plasma_complex64_t *A2, int lda2,
                       const plasma_complex64_t *V, int ldv,
                       const plasma_complex64_t *T, int ldt,
                       plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void core_ompss_ztsqrt(int m, int n, int ib,
                       plasma_complex64_t *A1, int lda1,
                       plasma_complex64_t *A2, int lda2,
                       plasma_complex64_t *T,  int ldt,
                       plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void core_ompss_zunmqr(plasma_enum_t side, plasma_enum_t trans,
                       int m, int n, int k, int ib,
                       const plasma_complex64_t *A, int lda,
                       const plasma_complex64_t *T, int ldt,
                             plasma_complex64_t *C, int ldc,
                       plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

#undef COMPLEX

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // ICL_CORE_BLAS_Z_H
