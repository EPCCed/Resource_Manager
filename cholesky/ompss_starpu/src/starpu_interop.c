#include <starpu.h>
#include <starpurm.h>
#include <morse.h>
#include "../include/interop.h"
#include "../include/ddss.h"

void interop_starpu_init()
{
	starpurm_initialize();
	starpurm_set_drs_enable(NULL);
	MORSE_Init(-1, -1);
	
	gen_cpu_masks(4, -1);
}

void interop_starpu_finalize()
{
	MORSE_Finalize();
	starpurm_shutdown();
}

static void wake_up_thread(void *args)
{
	nanos_unblock_task(args);
}

static void offload_starpu_dgemm(void *void_args)
{
	struct args_dgemm *args = (struct args_dgemm *)void_args;
	
	MORSE_dgemm(args->TRANS_A, args->TRANS_B, args->M, args->N, args->K,
						args->ALPHA, args->A, args->LDA, args->B, args->LDB,
						args->BETA, args->C, args->LDC);
}

void interop_dgemm( enum DDSS_TRANS TRANS_A, enum DDSS_TRANS TRANS_B,
	    int M, int N, int K,
        double ALPHA, double *A, int LDA,
                      double *B, int LDB,
        double BETA,  double *C, int LDC )
{
	struct args_dgemm *func_args = malloc(sizeof(*func_args));
	
	func_args->TRANS_A = TRANS_A;
	func_args->TRANS_B = TRANS_B;
	func_args->M = M;
	func_args->N = N;
	func_args->K = K;
	func_args->ALPHA = ALPHA;
	func_args->A = A;
	func_args->LDA = LDA;
	func_args->B = B;
	func_args->LDB = LDB;
	func_args->BETA = BETA;
	func_args->C = C;
	func_args->LDC = LDC;

	void *blocking_context = nanos_get_current_blocking_context();
	starpurm_spawn_kernel_on_cpus_callback(NULL, &offload_starpu_dgemm, func_args, get_hwloc_cpu_mask(0), wake_up_thread, blocking_context);
	nanos_block_current_task(blocking_context);
}

static void offload_starpu_dtrsm(void *void_args)
{
	struct args_dtrsm *args = (struct args_dtrsm *)void_args;
	
	MORSE_dtrsm(args->SIDE, args->UPLO, args->TRANS_A, args->DIAG,
						args->M, args->N, args->ALPHA,
						args->A, args->LDA,
						args->B, args->LDB);
}

void interop_dtrsm( enum DDSS_SIDE SIDE, enum DDSS_UPLO UPLO, 
		enum DDSS_TRANS TRANS_A,
        enum DDSS_DIAG DIAG, int M, int N,
        const double ALPHA, double *A, int LDA,
                            double *B, int LDB)
{
	struct args_dtrsm *func_args = malloc(sizeof(*func_args));
	
	func_args->SIDE = SIDE;
	func_args->UPLO = UPLO;
	func_args->TRANS_A = TRANS_A;
	func_args->DIAG = DIAG;
	func_args->M = M;
	func_args->N = N;
	func_args->ALPHA = ALPHA;
	func_args->A = A;
	func_args->LDA = LDA;
	func_args->B = B;
	func_args->LDB = LDB;

	void *blocking_context = nanos_get_current_blocking_context();
	starpurm_spawn_kernel_on_cpus_callback(NULL, &offload_starpu_dtrsm, func_args, get_hwloc_cpu_mask(1), wake_up_thread, blocking_context);
	nanos_block_current_task(blocking_context);
}

static void offload_starpu_dsyrk(void *void_args)
{
	struct args_dsyrk *args = (struct args_dsyrk *)void_args;
	
	MORSE_dsyrk(args->UPLO, args->TRANS, args->N, args->K,
						args->ALPHA, args->A, args->LDA,
						args->BETA, args->C, args->LDC);
}

void interop_dsyrk( enum DDSS_UPLO UPLO, enum DDSS_TRANS TRANS, 
		int N, int K,
        const double ALPHA, double *A, int LDA,
    	const double BETA,  double *C, int LDC)
{
	struct args_dsyrk *func_args = malloc(sizeof(*func_args));
	
	func_args->UPLO = UPLO;
	func_args->TRANS = TRANS;
	func_args->N = N;
	func_args->K = K;
	func_args->ALPHA = ALPHA;
	func_args->A = A;
	func_args->LDA = LDA;
	func_args->BETA = BETA;
	func_args->C = C;
	func_args->LDC = LDC;

	void *blocking_context = nanos_get_current_blocking_context();
	starpurm_spawn_kernel_on_cpus_callback(NULL, &offload_starpu_dsyrk, func_args, get_hwloc_cpu_mask(2), wake_up_thread, blocking_context);
	nanos_block_current_task(blocking_context);
}

static void offload_starpu_dpotrf(void *void_args)
{
	struct args_dpotrf *args = (struct args_dpotrf *)void_args;
	
	MORSE_dpotrf(args->UPLO, args->N, args->A, args->LDA);
}

void interop_dpotrf( enum DDSS_UPLO UPLO, int N, double *A, int LDA )
{
	struct args_dpotrf *func_args = malloc(sizeof(*func_args));
	
	func_args->UPLO = UPLO;
	func_args->N = N;
	func_args->A = A;
	func_args->LDA = LDA;

	void *blocking_context = nanos_get_current_blocking_context();
	starpurm_spawn_kernel_on_cpus_callback(NULL, &offload_starpu_dpotrf, func_args, get_hwloc_cpu_mask(3), wake_up_thread, blocking_context);
	nanos_block_current_task(blocking_context);
	
}
