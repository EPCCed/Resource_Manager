/**
 *
 * @fil
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
 * @ingroup core_geqrt
 *
 *  Computes a QR factorization of an m-by-n tile A:
 *  The factorization has the form
 *    \f[
 *        A = Q \times R
 *    \f]
 *  The tile Q is represented as a product of elementary reflectors
 *    \f[
 *        Q = H(1) H(2) ... H(k),
 *    \f]
 *  where \f$ k = min(m,n) \f$.
 *
 *  Each \f$ H(i) \f$ has the form
 *    \f[
 *        H(i) = I - \tau \times v \times v^H
 *    \f]
 *  where \f$ tau \f$ is a scalar, and \f$ v \f$ is a vector with
 *  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
 *  and \f$ tau \f$ in tau(i).
 *
 *******************************************************************************
 *
 * @param[in] m
 *          The number of rows of the tile A.  m >= 0.
 *
 * @param[in] n
 *         The number of columns of the tile A.  n >= 0.
 *
 * @param[in] ib
 *         The inner-blocking size.  ib >= 0.
 *
 * @param[in,out] A
 *         On entry, the m-by-n tile A.
 *         On exit, the elements on and above the diagonal of the array
 *         contain the min(m,n)-by-n upper trapezoidal tile R (R is
 *         upper triangular if m >= n); the elements below the diagonal,
 *         with the array tau, represent the unitary tile Q as a
 *         product of elementary reflectors (see Further Details).
 *
 * @param[in] lda
 *         The leading dimension of the array A.  lda >= max(1,m).
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
 * @param[in] lwork
 *         Size of the array work. Should be at least ib*n.
 *
 *******************************************************************************
 *
 * @retval PlasmaSuccess successful exit
 * @retval < 0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
__attribute__((weak))
int core_zgeqrt(int m, int n, int ib,
                plasma_complex64_t *A, int lda,
                plasma_complex64_t *T, int ldt,
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
    if ((ib < 0) || ( (ib == 0) && (m > 0) && (n > 0) )) {
        coreblas_error("illegal value of ib");
        return -3;
    }
    if (A == NULL) {
        coreblas_error("NULL A");
        return -4;
    }
    if (lda < imax(1, m) && m > 0) {
        coreblas_error("illegal value of lda");
        return -5;
    }
    if (T == NULL) {
        coreblas_error("NULL T");
        return -6;
    }
    if (ldt < imax(1, ib) && ib > 0) {
        coreblas_error("illegal value of ldt");
        return -7;
    }
    if (tau == NULL) {
        coreblas_error("NULL tau");
        return -8;
    }
    if (work == NULL) {
        coreblas_error("NULL work");
        return -9;
    }

    // quick return
    if (m == 0 || n == 0 || ib == 0)
        return PlasmaSuccess;

    int k = imin(m, n);
    for (int i = 0; i < k; i += ib) {
        int sb = imin(ib, k-i);

        LAPACKE_zgeqr2_work(LAPACK_COL_MAJOR,
                            m-i, sb,
                            &A[lda*i+i], lda,
                            &tau[i], work);

        LAPACKE_zlarft_work(LAPACK_COL_MAJOR,
                            lapack_const(PlasmaForward),
                            lapack_const(PlasmaColumnwise),
                            m-i, sb,
                            &A[lda*i+i], lda,
                            &tau[i],
                            &T[ldt*i], ldt);

        if (n > i+sb) {
            LAPACKE_zlarfb_work(LAPACK_COL_MAJOR,
                                lapack_const(PlasmaLeft),
                                lapack_const(Plasma_ConjTrans),
                                lapack_const(PlasmaForward),
                                lapack_const(PlasmaColumnwise),
                                m-i, n-i-sb, sb,
                                &A[lda*i+i],      lda,
                                &T[ldt*i],        ldt,
                                &A[lda*(i+sb)+i], lda,
                                work, n-i-sb);
        }
    }

    return PlasmaSuccess;
}

/******************************************************************************/
// The function to be run as a task.
static void core_starpu_cpu_zgeqrt(void *descr[], void *cl_arg)
{
    int m, n, ib;
    plasma_complex64_t *A, *T; 
    plasma_workspace_t work;
    int lda, ldt;

    // Unpack data of tiles.
    A = (plasma_complex64_t *) STARPU_MATRIX_GET_PTR(descr[0]);
    T = (plasma_complex64_t *) STARPU_MATRIX_GET_PTR(descr[1]);

    // Unpack scalar parameters.
    starpu_codelet_unpack_args(cl_arg, &m, &n, &ib, &lda, &ldt, &work);

    // Prepare workspaces.
    int id = starpu_worker_get_id();
    plasma_complex64_t *tau = (plasma_complex64_t*) work.spaces[id];
    
    // Call the kernel.
    core_zgeqrt(m, n, ib,
                A, lda,
                T, ldt,
                tau,
                tau+n); // work for LAPACK
}

/******************************************************************************/
// StarPU codelet.
struct starpu_codelet core_starpu_codelet_zgeqrt = {
    .cpu_func  = core_starpu_cpu_zgeqrt,
    .nbuffers  = 2,
    .name      = "zgeqrt"
};

/******************************************************************************/
// The function for inserting a task.
void core_starpu_zgeqrt(
    int m, int n, int ib,
    starpu_data_handle_t A, int lda,
    starpu_data_handle_t T, int ldt,
    plasma_workspace_t work, 
    plasma_sequence_t *sequence, plasma_request_t *request)
{
    starpu_insert_task(
        &core_starpu_codelet_zgeqrt,
        STARPU_VALUE,    &m,                 sizeof(int),
        STARPU_VALUE,    &n,                 sizeof(int),
        STARPU_VALUE,    &ib,                sizeof(int),
        STARPU_RW,       A,
        STARPU_VALUE,    &lda,               sizeof(int),
        STARPU_W,        T,
        STARPU_VALUE,    &ldt,               sizeof(int),
        STARPU_VALUE,    &work,              sizeof(plasma_workspace_t),
        STARPU_NAME, "zgeqrt",
        0);
}
