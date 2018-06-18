#ifndef STARPU_DGEMM_H
#define STARPU_DGEMM_H

void interop_starpu_init();

int interop_starpu_dgemm(int transa, int transb, int ydim, int xdim, int zdim,
			double alpha, double *A, int lda, double *B, int ldb, 
			double beta, double *C, int ldc);

void interop_starpu_finalize();

#endif
