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
#include "plasma_descriptor.h"
#include "plasma_internal.h"
#include "plasma_types.h"
#include "core_blas.h"

#define A(m, n) (plasma_complex64_t*)plasma_tile_addr(A, m, n)

/******************************************************************************/
void plasma_pzgeswp(plasma_enum_t colrow,
                    plasma_desc_t A, int *ipiv, int incx,
                    plasma_sequence_t *sequence, plasma_request_t *request)
{
    // Return if failed sequence.
    if (sequence->status != PlasmaSuccess)
        return;

    if (colrow == PlasmaRowwise) {
        for (int n = 0; n < A.nt; n++) {
            plasma_complex64_t *a00, *a10;

            a00 = A(0, n);
            a10 = A(A.mt-1, n);

            // Multidependency of the whole panel on its individual tiles.
            for (int m = 1; m < A.mt-1; m++) {
                plasma_complex64_t *amn = A(m, n);
                #pragma oss task in(amn[0]) \
                                 inout(a00[0]) \
                                 no_copy_deps
                {
                    int l = 1;
                    l++;
                }
            }

            int ma00 = (A.mt-1)*A.mb;
            int na00 = plasma_tile_nmain(A, n);

            int lda10 = plasma_tile_mmain(A, A.mt-1);
            int nva10 = plasma_tile_nview(A, n);

            #pragma oss task in(ipiv[0;A.m]) \
                             inout(a00[0;ma00*na00]) \
                             inout(a10[0;lda10*nva10]) \
                             no_copy_deps
            {
                int nvan = plasma_tile_nview(A, n);
                plasma_desc_t view = plasma_desc_view(A, 0, n*A.nb, A.m, nvan);
                core_zgeswp(colrow, view, 1, A.m, ipiv, incx);
            }

            // Multidependency of individual tiles on the whole panel.
            for (int m = 1; m < A.mt-1; m++) {
                plasma_complex64_t *amn = A(m, n);
                #pragma oss task in(a00[0]) \
                                 inout(amn[0]) \
                                 no_copy_deps
                {
                    int l = 1;
                    l++;
                }
            }
        }
    }
    else { // PlasmaColumnwise
        for (int m = 0; m < A.mt; m++) {
            plasma_complex64_t *a00, *a01;

            a00 = A(m, 0);
            a01 = A(m, A.nt-1);

            // Multidependency of the whole (row) panel on its individual tiles.
            for (int n = 1; n < A.nt-1; n++) {
                plasma_complex64_t *amn = A(m, n);
                #pragma oss task in(amn[0]) \
                                 inout(a00[0]) \
                                 no_copy_deps
                {
                    int l = 1;
                    l++;
                }
            }

            #pragma oss task in(ipiv[0;A.n]) \
                             inout(a00[0]) \
                             inout(a01[0]) \
                             no_copy_deps
            {
                int mvam = plasma_tile_mview(A, m);
                plasma_desc_t view = plasma_desc_view(A, m*A.mb, 0, mvam, A.n);
                core_zgeswp(colrow, view, 1, A.n, ipiv, incx);
            }

            // Multidependency of individual tiles on the whole (row) panel.
            for (int n = 1; n < A.nt-1; n++) {
                plasma_complex64_t *amn = A(m, n);
                #pragma oss task in(a00[0]) \
                                 inout(amn[0]) \
                                 no_copy_deps
                {
                    int l = 1;
                    l++;
                }
            }
        }
    }
}
