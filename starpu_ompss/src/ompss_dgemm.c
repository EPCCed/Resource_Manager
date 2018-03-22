#define _GNU_SOURCE
#include <sched.h>
#include <pthread.h>

#include <nanos6.h>
#include <nanos6/library-mode.h>

#include <ddss.h>

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
	int res;
};

struct ompss_cond_var {
	pthread_mutex_t mutex;
	pthread_cond_t cond;
	int signaled;
};

static void offload_ompss_dgemm(void *void_args)
{
	struct args_dgemm *args = (struct args_dgemm *)void_args;
	
	args->res = ddss_dgemm(args->TRANS_A, args->TRANS_B, args->M, args->N, args->K,
						args->ALPHA, args->A, args->LDA, args->B, args->LDB,
						args->BETA, args->C, args->LDC);
}

static void nanos6_wait_callback(void *args) {
	struct ompss_cond_var *cond_var = (struct ompss_cond_var *)args;
	pthread_mutex_lock(&cond_var->mutex);
	cond_var->signaled = 1;
	pthread_cond_signal(&cond_var->cond);
	pthread_mutex_unlock(&cond_var->mutex);
}

int interop_ompss_init(int ncpu)
{
	nanos_preinit();
	nanos_init();
	
	return 0;
}

int interop_ompss_dgemm(cpu_set_t *cpu_mask,
		int TRANS_A, int TRANS_B,
	    int M, int N, int K,
        double ALPHA, double *A, int LDA,
                      double *B, int LDB,
        double BETA,  double *C, int LDC)
{
	struct ompss_cond_var cond_var = {PTHREAD_MUTEX_INITIALIZER, PTHREAD_COND_INITIALIZER, 0};
	
	struct args_dgemm args;
	args.TRANS_A = TRANS_A;
	args.TRANS_B = TRANS_B;
	args.M = M;
	args.N = N;
	args.K = K;
	args.ALPHA = ALPHA;
	args.A = A;
	args.LDA = LDA;
	args.B = B;
	args.LDB = LDB;
	args.BETA = BETA;
	args.C = C;
	args.LDC = LDC;
	args.res = -1;

	nanos_spawn_function(&offload_ompss_dgemm, &args, &nanos6_wait_callback, &cond_var, "offload", cpu_mask);
	
	pthread_mutex_lock(&cond_var.mutex);
	while (cond_var.signaled == 0) {
		pthread_cond_wait(&cond_var.cond, &cond_var.mutex);
	}
	pthread_mutex_unlock(&cond_var.mutex);
	
	return args.res;
}

int interop_ompss_finalize()
{
	nanos_shutdown();
	
	return 0;
}
