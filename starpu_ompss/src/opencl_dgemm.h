#ifndef OPENCL_DGEMM_H
#define OPENCL_DGEMM_H

#include <CL/cl.h>

int interop_opencl_init(int max_cpus_task);

int interop_opencl_dgemm(int ctx_id,
		int TRANS_A, int TRANS_B,
	    int M, int N, int K,
        double ALPHA, double *A, int LDA,
                      double *B, int LDB,
        double BETA,  double *C, int LDC );

void interop_opencl_finalize();

#endif
