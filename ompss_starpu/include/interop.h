#ifndef INTEROP_H
#define INTEROP_H

#include <hwloc.h>
#include "ddss.h"

struct args_dgemm
{
	int TRANS_A;
	int TRANS_B;
	int M;
	int N;
	int K;
	double ALPHA;
	double *A;
	int LDA;
	double *B;
	int LDB;
	double BETA;
	double *C;
	int LDC;
};

struct args_dtrsm
{
	int SIDE;
	int UPLO;
	int TRANS_A;
	int DIAG;
	int M;
	int N;
	double ALPHA;
	double *A;
	int LDA;
	double *B;
	int LDB;
};

struct args_dsyrk
{
	int UPLO;
	int TRANS;
	int N;
	int K;
	double ALPHA;
	double *A;
	int LDA;
	double BETA;
	double *C;
	int LDC;
};

struct args_dpotrf
{
	int UPLO;
	int N;
	double *A;
	int LDA;
};

void interop_starpu_init();
void interop_starpu_finalize();

void interop_dgemm( enum DDSS_TRANS TRANS_A, enum DDSS_TRANS TRANS_B,
	    int M, int N, int K,
        double ALPHA, double *A, int LDA,
                      double *B, int LDB,
        double BETA,  double *C, int LDC );

void interop_dtrsm( enum DDSS_SIDE SIDE, enum DDSS_UPLO UPLO, 
		enum DDSS_TRANS TRANS_A,
        enum DDSS_DIAG DIAG, int M, int N,
        const double ALPHA, double *A, int LDA,
                            double *B, int LDB);

void interop_dsyrk( enum DDSS_UPLO UPLO, enum DDSS_TRANS TRANS, 
		int N, int K,
        const double ALPHA, double *A, int LDA,
    	const double BETA,  double *C, int LDC);

void interop_dpotrf( enum DDSS_UPLO UPLO, int N, double *A, int LDA );

void gen_cpu_masks(int num_kernels, int max_cpus_per_kernel);
cpu_set_t *get_cpu_mask(int kernel);
hwloc_cpuset_t get_hwloc_cpu_mask(int kernel);

#endif
