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

// This will be swapped during the automatic code generation.
#undef REAL
#define COMPLEX

/***************************************************************************//**
 *
 * @ingroup core_tsqrt
 *
 * Computes a QR factorization of a rectangular matrix
 * formed by coupling an n-by-n upper triangular tile A1
 * on top of an m-by-n tile A2:
 *
 *    | A1 | = Q * R
 *    | A2 |
 *
 *******************************************************************************
 *
 * @param[in] m
 *         The number of columns of the tile A2. m >= 0.
 *
 * @param[in] n
 *         The number of rows of the tile A1.
 *         The number of columns of the tiles A1 and A2. n >= 0.
 *
 * @param[in] ib
 *         The inner-blocking size.  ib >= 0.
 *
 * @param[in,out] A1
 *         On entry, the n-by-n tile A1.
 *         On exit, the elements on and above the diagonal of the array
 *         contain the n-by-n upper trapezoidal tile R;
 *         the elements below the diagonal are not referenced.
 *
 * @param[in] lda1
 *         The leading dimension of the array A1. LDA1 >= max(1,N).
 *
 * @param[in,out] A2
 *         On entry, the m-by-n tile A2.
 *         On exit, all the elements with the array tau, represent
 *         the unitary tile Q as a product of elementary reflectors
 *         (see Further Details).
 *
 * @param[in] lda2
 *         The leading dimension of the tile A2. lda2 >= max(1,m).
 *
 * @param[out] T
 *         The ib-by-n triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] ldt
 *         The leading dimension of the array T. ldt >= ib.
 *
 * @param tau
 *         Auxiliary workspace array of length n.
 *
 * @param work
 *         Auxiliary workspace array of length ib*n.
 *
 *******************************************************************************
 *
 * @retval PlasmaSuccess successful exit
 * @retval < 0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
__attribute__((weak))
int core_ztsqrt(int m, int n, int ib,
                plasma_complex64_t *A1, int lda1,
                plasma_complex64_t *A2, int lda2,
                plasma_complex64_t *T,  int ldt,
                plasma_complex64_t *tau,
                plasma_complex64_t *work)
{
    // Check input arguments.
    if (m < 0) {
        coreblas_error("illegal value of m");
        return -1;
    }
    if (n < 0) {
        coreblas_error("illegal value of n");
        return -2;
    }
    if (ib < 0) {
        coreblas_error("illegal value of ib");
        return -3;
    }
    if (A1 == NULL) {
        coreblas_error("NULL A1");
        return -4;
    }
    if (lda1 < imax(1, m) && m > 0) {
        coreblas_error("illegal value of lda1");
        return -5;
    }
    if (A2 == NULL) {
        coreblas_error("NULL A2");
        return -6;
    }
    if (lda2 < imax(1, m) && m > 0) {
        coreblas_error("illegal value of lda2");
        return -7;
    }
    if (T == NULL) {
        coreblas_error("NULL T");
        return -8;
    }
    if (ldt < imax(1, ib) && ib > 0) {
        coreblas_error("illegal value of ldt");
        return -9;
    }
    if (tau == NULL) {
        coreblas_error("NULL tau");
        return -10;
    }
    if (work == NULL) {
        coreblas_error("NULL work");
        return -11;
    }

    // quick return
    if (m == 0 || n == 0 || ib == 0)
        return PlasmaSuccess;

    static plasma_complex64_t zone  = 1.0;
    static plasma_complex64_t zzero = 0.0;

    for (int ii = 0; ii < n; ii += ib) {
        int sb = imin(n-ii, ib);
        for (int i = 0; i < sb; i++) {
            // Generate elementary reflector H( II*IB+I ) to annihilate
            // A( II*IB+I:M, II*IB+I ).
            LAPACKE_zlarfg_work(m+1, &A1[lda1*(ii+i)+ii+i], &A2[lda2*(ii+i)], 1,
                                &tau[ii+i]);

            if (ii+i+1 < n) {
                // Apply H( II*IB+I ) to A( II*IB+I:M, II*IB+I+1:II*IB+IB )
                // from the left.
                plasma_complex64_t alpha = -conj(tau[ii+i]);
                cblas_zcopy(sb-i-1, &A1[lda1*(ii+i+1)+(ii+i)], lda1, work, 1);
#ifdef COMPLEX
                LAPACKE_zlacgv_work(sb-i-1, work, 1);
#endif
                cblas_zgemv(CblasColMajor, (CBLAS_TRANSPOSE)Plasma_ConjTrans,
                            m, sb-i-1,
                            CBLAS_SADDR(zone), &A2[lda2*(ii+i+1)], lda2,
                            &A2[lda2*(ii+i)], 1,
                            CBLAS_SADDR(zone), work, 1);
#ifdef COMPLEX
                LAPACKE_zlacgv_work(sb-i-1, work, 1);
#endif
                cblas_zaxpy(sb-i-1, CBLAS_SADDR(alpha), work, 1,
                            &A1[lda1*(ii+i+1)+ii+i], lda1);
#ifdef COMPLEX
                LAPACKE_zlacgv_work(sb-i-1, work, 1);
#endif
                cblas_zgerc(CblasColMajor,
                            m, sb-i-1,
                            CBLAS_SADDR(alpha), &A2[lda2*(ii+i)], 1,
                            work, 1,
                            &A2[lda2*(ii+i+1)], lda2);
            }
            // Calculate T.
            plasma_complex64_t alpha = -tau[ii+i];
            cblas_zgemv(CblasColMajor, (CBLAS_TRANSPOSE)Plasma_ConjTrans,
                        m, i,
                        CBLAS_SADDR(alpha), &A2[lda2*ii], lda2,
                        &A2[lda2*(ii+i)], 1,
                        CBLAS_SADDR(zzero), &T[ldt*(ii+i)], 1);

            cblas_ztrmv(CblasColMajor, (CBLAS_UPLO)PlasmaUpper,
                        (CBLAS_TRANSPOSE)PlasmaNoTrans,
                        (CBLAS_DIAG)PlasmaNonUnit,
                        i,
                        &T[ldt*ii], ldt,
                        &T[ldt*(ii+i)], 1);

            T[ldt*(ii+i)+i] = tau[ii+i];
        }
        if (n > ii+sb) {
            core_ztsmqr(PlasmaLeft, Plasma_ConjTrans,
                        sb, n-(ii+sb), m, n-(ii+sb), ib, ib,
                        &A1[lda1*(ii+sb)+ii], lda1,
                        &A2[lda2*(ii+sb)], lda2,
                        &A2[lda2*ii], lda2,
                        &T[ldt*ii], ldt,
                        work, sb);
        }
    }

    return PlasmaSuccess;
}

/******************************************************************************/
// The function to be run as a task.
static void core_starpu_cpu_ztsqrt(void *descr[], void *cl_arg)
{
    int m, n, ib;
    plasma_complex64_t *A1, *A2, *T;
    int lda1, lda2, ldt;
    plasma_workspace_t work;

    // Unpack data of tiles.
    A1 = (plasma_complex64_t *) STARPU_MATRIX_GET_PTR(descr[0]);
    A2 = (plasma_complex64_t *) STARPU_MATRIX_GET_PTR(descr[1]);
    T  = (plasma_complex64_t *) STARPU_MATRIX_GET_PTR(descr[2]);

    // Unpack scalar parameters.
    starpu_codelet_unpack_args(cl_arg, &m, &n, &ib,
                               &lda1, &lda2, &ldt, &work);

    // Prepare workspaces.
    int id = starpu_worker_get_id();
    plasma_complex64_t *tau = (plasma_complex64_t*) work.spaces[id];

    // Call the kernel.
    core_ztsqrt(m, n, ib,
                A1, lda1,
                A2, lda2,
                T,  ldt,
                tau,
                tau+n);
}

/******************************************************************************/
// StarPU codelet.
struct starpu_codelet core_starpu_codelet_ztsqrt = {
    .cpu_func  = core_starpu_cpu_ztsqrt,
    .nbuffers  = 3,
    .name      = "ztsqrt"
};

/******************************************************************************/
void core_starpu_ztsqrt(int m, int n, int ib,
                        starpu_data_handle_t A1, int lda1,
                        starpu_data_handle_t A2, int lda2,
                        starpu_data_handle_t T,  int ldt,
                        plasma_workspace_t work,
                        plasma_sequence_t *sequence, plasma_request_t *request)
{
    starpu_insert_task(
        &core_starpu_codelet_ztsqrt,
        STARPU_VALUE,    &m,                 sizeof(int),
        STARPU_VALUE,    &n,                 sizeof(int),
        STARPU_VALUE,    &ib,                sizeof(int),
        STARPU_RW,       A1,
        STARPU_VALUE,    &lda1,              sizeof(int),
        STARPU_RW,       A2,
        STARPU_VALUE,    &lda2,              sizeof(int),
        STARPU_W,        T,
        STARPU_VALUE,    &ldt,               sizeof(int),
        STARPU_VALUE,    &work,              sizeof(plasma_workspace_t),
        STARPU_NAME, "ztsqrt",
        0);
}
