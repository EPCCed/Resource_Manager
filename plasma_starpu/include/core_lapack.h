/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 **/
#ifndef ICL_CORE_LAPACK_H
#define ICL_CORE_LAPACK_H

#if defined(HAVE_MKL) || defined(PLASMA_WITH_MKL)
    #define MKL_Complex16 double _Complex
    #define MKL_Complex8  float _Complex

    #include <mkl_cblas.h>
    #include <mkl_lapacke.h>

    // MKL LAPACKE doesn't provide LAPACK_GLOBAL macro, so define it here.
    // MKL provides all 3 name manglings (foo, foo_, FOO); pick foo_.
    #ifndef LAPACK_GLOBAL
    #define LAPACK_GLOBAL(lcname,UCNAME)  lcname##_
    #endif
#elif defined(HAVE_ESSL) || defined(PLASMA_WITH_ESSL)
    // GCC + ESSL(BLAS) + LAPACKE/CBLAS from LAPACK
    #include <cblas.h>
    #include <lapacke.h>

    #ifndef LAPACK_GLOBAL
    #define LAPACK_GLOBAL(lcname,UCNAME) lcname##_
    #endif
#else
    #include <cblas.h>
    #include <lapacke.h>

    // Intel mkl_cblas.h does: typedef enum {...} CBLAS_ORDER;
    // Netlib    cblas.h does: enum CBLAS_ORDER {...};
    // OpenBLAS  cblas.h does: typedef enum CBLAS_ORDER {...} CBLAS_ORDER;
    // We use (CBLAS_ORDER), so add these typedefs for Netlib.
    #ifndef OPENBLAS_VERSION
    typedef enum CBLAS_ORDER CBLAS_ORDER;
    typedef enum CBLAS_TRANSPOSE CBLAS_TRANSPOSE;
    typedef enum CBLAS_UPLO CBLAS_UPLO;
    typedef enum CBLAS_DIAG CBLAS_DIAG;
    typedef enum CBLAS_SIDE CBLAS_SIDE;
    #endif
#endif

#include "core_lapack_s.h"
#include "core_lapack_d.h"
#include "core_lapack_c.h"
#include "core_lapack_z.h"

#endif // ICL_CORE_LAPACK_H
