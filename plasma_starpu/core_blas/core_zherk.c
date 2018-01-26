/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @precisions normal z -> c
 *
 **/

#include "core_blas.h"
#include "plasma_types.h"
#include "core_lapack.h"

#include "starpu.h"

/***************************************************************************//**
 *
 * @ingroup core_herk
 *
 *  Performs one of the Hermitian rank k operations
 *
 *    \f[ C = \alpha A \times A^H + \beta C, \f]
 *    or
 *    \f[ C = \alpha A^H \times A + \beta C, \f]
 *
 *  where alpha and beta are real scalars, C is an n-by-n Hermitian
 *  matrix, and A is an n-by-k matrix in the first case and a k-by-n
 *  matrix in the second case.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          - PlasmaUpper: Upper triangle of C is stored;
 *          - PlasmaLower: Lower triangle of C is stored.
 *
 * @param[in] trans
 *          - PlasmaNoTrans:   \f[ C = \alpha A \times A^H + \beta C; \f]
 *          - PlasmaConjTrans: \f[ C = \alpha A^H \times A + \beta C. \f]
 *
 * @param[in] n
 *          The order of the matrix C. n >= 0.
 *
 * @param[in] k
 *          If trans = PlasmaNoTrans, number of columns of the A matrix;
 *          if trans = PlasmaConjTrans, number of rows of the A matrix.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 *          A is an lda-by-ka matrix.
 *          If trans = PlasmaNoTrans,   ka = k;
 *          if trans = PlasmaConjTrans, ka = n.
 *
 * @param[in] lda
 *          The leading dimension of the array A.
 *          If trans = PlasmaNoTrans,   lda >= max(1, n);
 *          if trans = PlasmaConjTrans, lda >= max(1, k).
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in,out] C
 *          C is an ldc-by-n matrix.
 *          On exit, the uplo part of the matrix is overwritten
 *          by the uplo part of the updated matrix.
 *
 * @param[in] ldc
 *          The leading dimension of the array C. ldc >= max(1, n).
 *
 ******************************************************************************/
__attribute__((weak))
void core_zherk(plasma_enum_t uplo, plasma_enum_t trans,
                int n, int k,
                double alpha, const plasma_complex64_t *A, int lda,
                double beta,        plasma_complex64_t *C, int ldc)
{
    cblas_zherk(CblasColMajor,
                (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                n, k,
                alpha, A, lda,
                beta,  C, ldc);
}

/******************************************************************************/
// StarPU HERK CPU kernel.
static void core_starpu_cpu_zherk(void *descr[], void * cl_arg)
{
    plasma_enum_t uplo,  trans;
    int n,  k;
    double alpha, beta;
    plasma_complex64_t *A;
    plasma_complex64_t *C;
    int lda, ldc;

    // Unpack data of tiles.
    A = (plasma_complex64_t *) STARPU_MATRIX_GET_PTR(descr[0]);
    C = (plasma_complex64_t *) STARPU_MATRIX_GET_PTR(descr[1]);

    // Unpack scalar parameters.
    starpu_codelet_unpack_args(cl_arg, &uplo, &trans, &n, &k, &alpha, &lda,
                               &beta, &ldc);

    // Call the kernel.
    core_zherk(uplo, trans,
               n, k,
               alpha, A, lda,
               beta,  C, ldc);
}

/******************************************************************************/
// StarPU codelet.
struct starpu_codelet core_starpu_codelet_zherk = {
    .cpu_func = core_starpu_cpu_zherk,
    .nbuffers = 2,
    .name     = "zherk"
};
    
/******************************************************************************/
// The function for task insertion.
void core_starpu_zherk(
    plasma_enum_t uplo, plasma_enum_t trans,
    int n, int k,
    double alpha, starpu_data_handle_t A, int lda,
    double beta,  starpu_data_handle_t C, int ldc,
    plasma_sequence_t *sequence, plasma_request_t *request)
{
    starpu_insert_task(
        &core_starpu_codelet_zherk,
        STARPU_VALUE,    &uplo,        sizeof(plasma_enum_t),
        STARPU_VALUE,    &trans,       sizeof(plasma_enum_t),
        STARPU_VALUE,    &n,           sizeof(int),
        STARPU_VALUE,    &k,           sizeof(int),
        STARPU_VALUE,    &alpha,       sizeof(double),
        STARPU_R,        A,
        STARPU_VALUE,    &lda,         sizeof(int),
        STARPU_VALUE,    &beta,        sizeof(double),
        STARPU_RW,       C,
        STARPU_VALUE,    &ldc,         sizeof(int),
        STARPU_NAME,     "zherk",
        0);
}
