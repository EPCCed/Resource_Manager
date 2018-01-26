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

#include "plasma_async.h"
#include "plasma_context.h"
#include "plasma_descriptor.h"
#include "plasma_internal.h"
#include "plasma_types.h"
#include "plasma_workspace.h"
#include "core_blas.h"

#define A(m, n) (plasma_complex64_t*)plasma_tile_addr(A, m, n)

/******************************************************************************/
void plasma_pzgetrf(plasma_desc_t A, int *ipiv,
                    plasma_sequence_t *sequence, plasma_request_t *request)
{
    // Return if failed sequence.
    if (sequence->status != PlasmaSuccess)
        return;

    // Read parameters from the context.
    plasma_context_t *plasma = plasma_context_self();

    // Set tiling parameters.
    int ib = plasma->ib;

    int minmtnt = imin(A.mt, A.nt);

    for (int k = 0; k < minmtnt; k++) {

        plasma_complex64_t *a00, *a20;
        a00 = A(k, k);
        a20 = A(A.mt-1, k);

        // Create fake dependencies of the whole panel on its individual tiles.
        // These tasks are inserted to generate a correct DAG rather than
        // doing any useful work.
        for (int m = k+1; m < A.mt-1; m++) {
            plasma_complex64_t *amk = A(m, k);
            #pragma oss task in(amk[0]) \
                             inout(a00[0]) \
                             priority(1) \
                             no_copy_deps
            {
                // Do some funny work here. It appears so that the compiler
                // might not insert the task if it is completely empty.
                int l = 1;
                l++;
            }
        }

        int ma00k = (A.mt-k-1)*A.mb;
        int na00k = plasma_tile_nmain(A, k);
        int lda20 = plasma_tile_mmain(A, A.mt-1);

        int nvak = plasma_tile_nview(A, k);
        int mvak = plasma_tile_mview(A, k);
        int ldak = plasma_tile_mmain(A, k);

        int num_panel_threads = imin(plasma->max_panel_threads,
                                     minmtnt-k);
        // panel
        #pragma oss task inout(a00[0;ma00k*na00k]) \
                         inout(a20[0;lda20*nvak]) \
                         out(ipiv[k*A.mb;mvak]) \
                         priority(1) \
                         no_copy_deps
        {
            volatile int *max_idx = (int*)malloc(num_panel_threads*sizeof(int));
            if (max_idx == NULL)
                plasma_request_fail(sequence, request, PlasmaErrorOutOfMemory);

            volatile plasma_complex64_t *max_val =
                (plasma_complex64_t*)malloc(num_panel_threads*sizeof(
                                            plasma_complex64_t));
            if (max_val == NULL)
                plasma_request_fail(sequence, request, PlasmaErrorOutOfMemory);

            volatile int info = 0;

            plasma_barrier_t barrier;
            plasma_barrier_init(&barrier);

            if (sequence->status == PlasmaSuccess) {
                // If nesting would not be expensive on architectures such as
                // KNL, this would resolve the issue with deadlocks caused by 
                // tasks expected to run are in fact not launched.
                //#pragma oss parallel for shared(barrier) 
                //                         schedule(dynamic,1) 
                //                         num_threads(num_panel_threads)
                //#pragma oss taskloop untied shared(barrier) 
                //                     num_tasks(num_panel_threads) 
                //                     priority(2)
                for (int rank = 0; rank < num_panel_threads; rank++) {
                    #pragma oss task shared(barrier) \
                                     priority(2) \
                                     no_copy_deps
                    {
                        plasma_desc_t view =
                            plasma_desc_view(A,
                                             k*A.mb, k*A.nb,
                                             A.m-k*A.mb, nvak);

                        core_zgetrf(view, &ipiv[k*A.mb], ib,
                                    rank, num_panel_threads,
                                    max_idx, max_val, &info,
                                    &barrier);

                        if (info != 0)
                            plasma_request_fail(sequence, request, k*A.mb+info);
                    }
                }
            }
            #pragma oss taskwait

            free((void*)max_idx);
            free((void*)max_val);

            for (int i = k*A.mb+1; i <= imin(A.m, k*A.mb+nvak); i++)
                ipiv[i-1] += k*A.mb;
        }
        
        // update
        for (int n = k+1; n < A.nt; n++) {
            plasma_complex64_t *a01, *a11, *a21;
            a01 = A(k, n);
            a11 = A(k+1, n);
            a21 = A(A.mt-1, n);

            int ma11k = (A.mt-k-2)*A.mb;
            int na11n = plasma_tile_nmain(A, n);
            int lda21 = plasma_tile_mmain(A, A.mt-1);

            int nvan = plasma_tile_nview(A, n);

            #pragma oss task in(a00[0;ma00k*na00k]) \
                             in(a20[0;lda20*nvak]) \
                             in(ipiv[k*A.mb;mvak]) \
                             inout(a01[0;ldak*nvan]) \
                             inout(a11[0;ma11k*na11n]) \
                             inout(a21[0;lda21*nvan]) \
                             priority(n == k+1) \
                             no_copy_deps
            {
                if (sequence->status == PlasmaSuccess) {
                    // geswp
                    int k1 = k*A.mb+1;
                    int k2 = imin(k*A.mb+A.mb, A.m);
                    plasma_desc_t view =
                        plasma_desc_view(A, 0, n*A.nb, A.m, nvan);
                    core_zgeswp(PlasmaRowwise, view, k1, k2, ipiv, 1);

                    // trsm
                    core_ztrsm(PlasmaLeft, PlasmaLower,
                               PlasmaNoTrans, PlasmaUnit,
                               mvak, nvan,
                               1.0, A(k, k), ldak,
                                    A(k, n), ldak);
                    // gemm
                    for (int m = k+1; m < A.mt; m++) {
                        int mvam = plasma_tile_mview(A, m);
                        int ldam = plasma_tile_mmain(A, m);

                        #pragma oss task priority(n == k+1) \
                                         no_copy_deps
                        {
                            core_zgemm(
                                PlasmaNoTrans, PlasmaNoTrans,
                                mvam, nvan, A.nb,
                                -1.0, A(m, k), ldam,
                                      A(k, n), ldak,
                                1.0,  A(m, n), ldam);
                        }
                    }
                }
                #pragma oss taskwait
            }
        }
    }

    // Multidependency of the whole ipiv on the individual chunks
    // corresponding to tiles. 
    for (int m = 0; m < minmtnt; m++) {
        // insert dummy task
        #pragma oss task in(ipiv[m*A.mb]) \
                         inout(ipiv[0]) \
                         no_copy_deps
        {
            int l = 1;
            l++;
        }
    }

    // pivoting to the left
    for (int k = 0; k < minmtnt-1; k++) {
        plasma_complex64_t *a10, *a20;
        a10 = A(k+1, k);
        a20 = A(A.mt-1, k);

        int ma10k = (A.mt-k-2)*A.mb;
        int na00k = plasma_tile_nmain(A, k);
        int lda20 = plasma_tile_mmain(A, A.mt-1);

        int nvak = plasma_tile_nview(A, k);

        #pragma oss task in(ipiv[0;imin(A.m,A.n)]) \
                         inout(a10[0;ma10k*na00k]) \
                         inout(a20[0;lda20*nvak]) \
                         no_copy_deps
        {
            if (sequence->status == PlasmaSuccess) {
                plasma_desc_t view =
                    plasma_desc_view(A, 0, k*A.nb, A.m, A.nb);
                int k1 = (k+1)*A.mb+1;
                int k2 = imin(A.m, A.n);
                core_zgeswp(PlasmaRowwise, view, k1, k2, ipiv, 1);
            }
        }

        // Multidependency of individual tiles on the whole panel.
        for (int m = k+2; m < A.mt-1; m++) {
            plasma_complex64_t *amk = A(m, k);
            #pragma oss task in(a10[0]) \
                             inout(amk[0]) \
                             no_copy_deps
            {
                // Do some funny work here. It appears so that the compiler
                // might not insert the task if it is completely empty.
                int l = 1;
                l++;
            }
        }
    }
}
