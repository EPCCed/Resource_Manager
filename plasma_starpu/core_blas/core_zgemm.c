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
#include "core_lapack.h"

#include "starpu.h"

/***************************************************************************//**
 *
 * @ingroup core_gemm
 *
 *  Performs one of the matrix-matrix operations
 *
 *    \f[ C = \alpha [op( A )\times op( B )] + \beta C, \f]
 *
 *  where op( X ) is one of:
 *    \f[ op( X ) = X,   \f]
 *    \f[ op( X ) = X^T, \f]
 *    \f[ op( X ) = X^H, \f]
 *
 *  alpha and beta are scalars, and A, B and C  are matrices, with op( A )
 *  an m-by-k matrix, op( B ) a k-by-n matrix and C an m-by-n matrix.
 *
 *******************************************************************************
 *
 * @param[in] transa
 *          - PlasmaNoTrans:   A is not transposed,
 *          - PlasmaTrans:     A is transposed,
 *          - PlasmaConjTrans: A is conjugate transposed.
 *
 * @param[in] transb
 *          - PlasmaNoTrans:   B is not transposed,
 *          - PlasmaTrans:     B is transposed,
 *          - PlasmaConjTrans: B is conjugate transposed.
 *
 * @param[in] m
 *          The number of rows of the matrix op( A ) and of the matrix C.
 *          m >= 0.
 *
 * @param[in] n
 *          The number of columns of the matrix op( B ) and of the matrix C.
 *          n >= 0.
 *
 * @param[in] k
 *          The number of columns of the matrix op( A ) and the number of rows
 *          of the matrix op( B ). k >= 0.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 *          An lda-by-ka matrix, where ka is k when transa = PlasmaNoTrans,
 *          and is m otherwise.
 *
 * @param[in] lda
 *          The leading dimension of the array A.
 *          When transa = PlasmaNoTrans, lda >= max(1,m),
 *          otherwise, lda >= max(1,k).
 *
 * @param[in] B
 *          An ldb-by-kb matrix, where kb is n when transb = PlasmaNoTrans,
 *          and is k otherwise.
 *
 * @param[in] ldb
 *          The leading dimension of the array B.
 *          When transb = PlasmaNoTrans, ldb >= max(1,k),
 *          otherwise, ldb >= max(1,n).
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in,out] C
 *          An ldc-by-n matrix. On exit, the array is overwritten by the m-by-n
 *          matrix ( alpha*op( A )*op( B ) + beta*C ).
 *
 * @param[in] ldc
 *          The leading dimension of the array C. ldc >= max(1,m).
 *
 ******************************************************************************/
__attribute__((weak))
void core_zgemm(plasma_enum_t transa, plasma_enum_t transb,
                int m, int n, int k,
                plasma_complex64_t alpha, const plasma_complex64_t *A, int lda,
                                          const plasma_complex64_t *B, int ldb,
                plasma_complex64_t beta,        plasma_complex64_t *C, int ldc)
{
    cblas_zgemm(CblasColMajor,
                (CBLAS_TRANSPOSE)transa, (CBLAS_TRANSPOSE)transb,
                m, n, k,
                CBLAS_SADDR(alpha), A, lda,
                                    B, ldb,
                CBLAS_SADDR(beta),  C, ldc);
}

/******************************************************************************/
// StarPU GEMM CPU kernel.
static void core_starpu_cpu_zgemm(void *descr[], void *cl_arg)
{
    plasma_enum_t transa, transb;
    int m, n, k;
    plasma_complex64_t alpha, beta;
    plasma_complex64_t *A, *B, *C;
    int lda, ldb, ldc;

    // Unpack data of tiles.
    A = (plasma_complex64_t *) STARPU_MATRIX_GET_PTR(descr[0]);
    B = (plasma_complex64_t *) STARPU_MATRIX_GET_PTR(descr[1]);
    C = (plasma_complex64_t *) STARPU_MATRIX_GET_PTR(descr[2]);

    // Unpack scalar parameters.
    starpu_codelet_unpack_args(cl_arg, &transa, &transb, &m, &n, &k, &alpha,
                               &lda, &ldb, &beta, &ldc);

    // Call the kernel.
    core_zgemm(transa, transb,
               m, n, k,
               alpha, A, lda,
               B, ldb,
               beta, C, ldc);
}

/******************************************************************************/
// StarPU codelet.
struct starpu_codelet core_starpu_codelet_zgemm = {
    .cpu_func  = core_starpu_cpu_zgemm,
    .nbuffers  = 3,
    .name      = "zgemm"
};       

/******************************************************************************/
// The function for task insertion.
void core_starpu_zgemm(
    plasma_enum_t transa, plasma_enum_t transb,
    int m, int n, int k,
    plasma_complex64_t alpha, starpu_data_handle_t A, int lda,
                              starpu_data_handle_t B, int ldb,
    plasma_complex64_t beta,  starpu_data_handle_t C, int ldc,
    plasma_sequence_t *sequence, plasma_request_t *request)
{
    starpu_insert_task(
        &core_starpu_codelet_zgemm,
        STARPU_VALUE,    &transa,            sizeof(plasma_enum_t),
        STARPU_VALUE,    &transb,            sizeof(plasma_enum_t),
        STARPU_VALUE,    &m,                 sizeof(int),
        STARPU_VALUE,    &n,                 sizeof(int),
        STARPU_VALUE,    &k,                 sizeof(int),
        STARPU_VALUE,    &alpha,             sizeof(plasma_complex64_t),
        STARPU_R,        A,
        STARPU_VALUE,    &lda,               sizeof(int),
        STARPU_R,        B,
        STARPU_VALUE,    &ldb,               sizeof(int),
        STARPU_VALUE,    &beta,              sizeof(plasma_complex64_t),
        STARPU_RW,       C,
        STARPU_VALUE,    &ldc,               sizeof(int),
        STARPU_NAME, "zgemm",
        0);
}
