#ifndef OMPSS_DGEMM_H
#define OMPSS_DGEMM_H

int interop_ompss_init();
int interop_ompss_dgemm(cpu_set_t *cpu_mask,
		int TRANS_A, int TRANS_B,
	    int M, int N, int K,
        double ALPHA, double *A, int LDA,
                      double *B, int LDB,
        double BETA,  double *C, int LDC);
int interop_ompss_finalize();

#endif
