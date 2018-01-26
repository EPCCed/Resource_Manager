/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @precisions normal z -> c d s
 *
 **/

#include "core_blas.h"
#include "plasma_types.h"
#include "plasma_internal.h"
#include "core_lapack.h"

#include "starpu.h"

/***************************************************************************//**
 *
 * @ingroup core_lacpy
 *
 *  Copies all or part of a two-dimensional matrix A to another matrix B.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          - PlasmaGeneral: entire A,
 *          - PlasmaUpper:   upper triangle,
 *          - PlasmaLower:   lower triangle.
 *
 * @param[in] transa
 *          - PlasmaNoTrans:   A is not transposed,
 *          - PlasmaTrans:     A is transposed,
 *          - PlasmaConjTrans: A is conjugate transposed.
 *
 * @param[in] m
 *          The number of rows of the matrices A and B.
 *          m >= 0.
 *
 * @param[in] n
 *          The number of columns of the matrices A and B.
 *          n >= 0.
 *
 * @param[in] A
 *          The m-by-n matrix to copy.
 *
 * @param[in] lda
 *          The leading dimension of the array A.
 *          lda >= max(1,m).
 *
 * @param[out] B
 *          The m-by-n copy of the matrix A.
 *          On exit, B = A ONLY in the locations specified by uplo.
 *
 * @param[in] ldb
 *          The leading dimension of the array B.
 *          ldb >= max(1,m).
 *
 ******************************************************************************/
__attribute__((weak))
void core_zlacpy(plasma_enum_t uplo, plasma_enum_t transa,
                 int m, int n,
                 const plasma_complex64_t *A, int lda,
                       plasma_complex64_t *B, int ldb)
{
    if (transa == PlasmaNoTrans) {
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR,
                            lapack_const(uplo),
                            m, n,
                            A, lda,
                            B, ldb);
    } else if (transa == PlasmaTrans) {
        switch(uplo) {
        case PlasmaUpper:
            for (int i = 0; i < imin(m, n); i++)
                for (int j = i; j < n; j++)
                    B[j + i*ldb] = A[i + j*lda];
            break;
        case PlasmaLower:
            for (int i = 0; i < m; i++)
                for (int j = 0; j <= imin(i, n); j++)
                    B[j + i*ldb] = A[i + j*lda];
            break;
        case PlasmaGeneral:
            for (int i = 0; i < m; i++)
                for (int j = 0; j < n; j++)
                    B[j + i*ldb] = A[i + j*lda];
            break;
        }
    } else {
        switch(uplo) {
        case PlasmaUpper:
            for (int i = 0; i < imin(m, n); i++)
                for (int j = 0; j < m; j++)
                    B[j + i*ldb] = conj(A[i + j*lda]);
            break;
        case PlasmaLower:
            for (int i = 0; i < m; i++)
                for (int j = 0; j <= imin(i, n); j++)
                    B[j + i*ldb] = conj(A[i + j*lda]);
            break;
        case PlasmaGeneral:
            for (int i = 0; i < m; i++)
                for (int j = 0; j < n; j++)
                    B[j + i*ldb] = conj(A[i + j*lda]);
            break;
        }
    }
}

/******************************************************************************/
// StarPU LACPY CPU kernel.
static void core_starpu_cpu_zlacpy_forward(void *descr[], void *cl_arg)
{
    plasma_enum_t uplo, transa;
    int x1, x2, y1, y2;
    plasma_complex64_t *A, *B;
    int lda, ldb;

    // Unpack data of tiles.
    B = (plasma_complex64_t *) STARPU_MATRIX_GET_PTR(descr[0]);

    // Unpack scalar parameters.
    starpu_codelet_unpack_args(cl_arg, &uplo, &transa,
                               &x1, &x2, &y1, &y2, 
                               &A, &lda, &ldb);

    // Call the kernel.
    core_zlacpy(uplo, transa,
                y2-y1, x2-x1,
                &(A[x1*lda+y1]), lda,
                &(B[x1*ldb+y1]), ldb);
}

static void core_starpu_cpu_zlacpy_backward(void *descr[], void *cl_arg)
{
    plasma_enum_t uplo, transa;
    int x1, x2, y1, y2;
    plasma_complex64_t *A, *B;
    int lda, ldb;

    // Unpack data of tiles.
    B = (plasma_complex64_t *) STARPU_MATRIX_GET_PTR(descr[0]);

    // Unpack scalar parameters.
    starpu_codelet_unpack_args(cl_arg, &uplo, &transa,
                               &x1, &x2, &y1, &y2, 
                               &A, &lda, &ldb);

    // Call the kernel.
    core_zlacpy(uplo, transa,
                y2-y1, x2-x1,
                &(B[x1*ldb+y1]), ldb,
                &(A[x1*lda+y1]), lda);
}

/******************************************************************************/
// StarPU codelets.
struct starpu_codelet core_starpu_codelet_zlacpy_forward = {
    .cpu_func  = core_starpu_cpu_zlacpy_forward,
    .nbuffers  = 1,
    .name      = "zlacpy_forward"
};       

struct starpu_codelet core_starpu_codelet_zlacpy_backward = {
    .cpu_func  = core_starpu_cpu_zlacpy_backward,
    .nbuffers  = 1,
    .name      = "zlacpy_backward"
};       

/******************************************************************************/
// The function for task insertion.
void core_starpu_zlacpy(
    plasma_enum_t uplo, plasma_enum_t transa, plasma_enum_t direction,
    int x1, int x2, int y1, int y2,
    plasma_complex64_t *A, int lda,
    starpu_data_handle_t B, int ldb,
    plasma_sequence_t *sequence, plasma_request_t *request)
{
    starpu_insert_task(
        (direction == PlasmaForward) ? 
            &core_starpu_codelet_zlacpy_forward : 
            &core_starpu_codelet_zlacpy_backward,
        STARPU_VALUE,    &uplo,              sizeof(plasma_enum_t),
        STARPU_VALUE,    &transa,            sizeof(plasma_enum_t),
        STARPU_VALUE,    &x1,                sizeof(int),
        STARPU_VALUE,    &x2,                sizeof(int),
        STARPU_VALUE,    &y1,                sizeof(int),
        STARPU_VALUE,    &y2,                sizeof(int),
        STARPU_VALUE,    &A,                 sizeof(plasma_complex64_t*),
        STARPU_VALUE,    &lda,               sizeof(int),
        (direction == PlasmaForward) ? STARPU_W : STARPU_R, B,
        STARPU_VALUE,    &ldb,               sizeof(int),
        STARPU_NAME, "zlacpy",
        0);
}
