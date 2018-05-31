/*
 * Based on StarPU matrix multiply example: starpu/examples/mult
 */
#define _GNU_SOURCE
#include <assert.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include <sched.h>
#include <sys/types.h>
#include <starpu.h>
#include <mkl.h>

#define KERNEL_BLAS 1
#define KERNEL_OMPSS 2
#define KERNEL_OPENCL 3

#ifdef ENABLE_OPENCL
#define INTEROP_KERNEL KERNEL_OPENCL
#else
#define INTEROP_KERNEL KERNEL_OMPSS
#endif

#if INTEROP_KERNEL == KERNEL_OMPSS
#include "ompss_dgemm.h"
#elif INTEROP_KERNEL == KERNEL_OPENCL
#include "opencl_dgemm.h"
#endif

#define STARPU_TILE_SIZE 4096
#define MAX_CPUS_TASK 4

struct starpu_dgemm_args
{
	double alpha;
	double beta;
#if INTEROP_KERNEL == KERNEL_OMPSS
	cpu_set_t cpu_mask;
#elif INTEROP_KERNEL == KERNEL_OPENCL
	int ctx_id;
#endif
};

static void partition_mult_data(int xdim, int xt, int ydim, int yt, int zdim,
								double *A, starpu_data_handle_t *A_handle, int lda,
								double *B, starpu_data_handle_t *B_handle, int ldb,
								double *C, starpu_data_handle_t *C_handle, int ldc)
{
	starpu_matrix_data_register(A_handle, STARPU_MAIN_RAM, (uintptr_t)A,
		lda, ydim, zdim, sizeof(double));
	starpu_matrix_data_register(B_handle, STARPU_MAIN_RAM, (uintptr_t)B,
		ldb, zdim, xdim, sizeof(double));
	starpu_matrix_data_register(C_handle, STARPU_MAIN_RAM, (uintptr_t)C,
		ldc, ydim, xdim, sizeof(double));

	struct starpu_data_filter vert;
	memset(&vert, 0, sizeof(vert));
	vert.filter_func = starpu_matrix_filter_vertical_block;
	vert.nchildren = xt;

	struct starpu_data_filter horiz;
	memset(&horiz, 0, sizeof(horiz));
	horiz.filter_func = starpu_matrix_filter_block;
	horiz.nchildren = yt;

	starpu_data_partition(*B_handle, &vert);
	starpu_data_partition(*A_handle, &horiz);

	starpu_data_map_filters(*C_handle, 2, &vert, &horiz);
}

void cpu_mult(void *descr[], void *arg)
{
	struct starpu_dgemm_args *params = (struct starpu_dgemm_args *)arg;
	
	double alpha = params->alpha;
	double beta = params->beta;
	
	double *subA = (double *)STARPU_MATRIX_GET_PTR(descr[0]);
	double *subB = (double *)STARPU_MATRIX_GET_PTR(descr[1]);
	double *subC = (double *)STARPU_MATRIX_GET_PTR(descr[2]);

	unsigned nxC = STARPU_MATRIX_GET_NX(descr[2]);
	unsigned nyC = STARPU_MATRIX_GET_NY(descr[2]);
	unsigned nyA = STARPU_MATRIX_GET_NY(descr[0]);

	unsigned ldA = STARPU_MATRIX_GET_LD(descr[0]);
	unsigned ldB = STARPU_MATRIX_GET_LD(descr[1]);
	unsigned ldC = STARPU_MATRIX_GET_LD(descr[2]);
	
#if INTEROP_KERNEL == KERNEL_BLAS
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
				nxC, nyC, nyA, alpha, subA, ldA, subB, ldB, beta, subC, ldC);
#elif INTEROP_KERNEL == KERNEL_OMPSS
	interop_ompss_dgemm(&params->cpu_mask, CblasNoTrans, CblasNoTrans,
				nxC, nyC, nyA, alpha, subA, ldA, subB, ldB, beta, subC, ldC);
#elif INTEROP_KERNEL == KERNEL_OPENCL
	interop_opencl_dgemm(params->ctx_id, CblasNoTrans, CblasNoTrans,
				nxC, nyC, nyA, alpha, subA, ldA, subB, ldB, beta, subC, ldC);
#endif
}

static struct starpu_codelet cl =
{
	.type = STARPU_SEQ, /* changed to STARPU_SPMD if -spmd is passed */
	.max_parallelism = INT_MAX,
	.cpu_funcs = {cpu_mult},
	.cpu_funcs_name = {"cpu_mult"},
	.nbuffers = 3,
	.modes = {STARPU_R, STARPU_R, STARPU_RW},
};

void interop_starpu_init(int num_kernels)
{
	int ret = starpu_init(NULL);
	assert(ret == 0);
}

int interop_starpu_dgemm(int transa, int transb, int ydim, int xdim, int zdim,
			double alpha, double *A, int lda, double *B, int ldb, 
			double beta, double *C, int ldc)
{
	int ret = 0;
	int xt, yt;
	starpu_data_handle_t A_handle, B_handle, C_handle;
	
	if (xdim % STARPU_TILE_SIZE == 0) {
		xt = xdim / STARPU_TILE_SIZE;
	} else {
		xt = (xdim / STARPU_TILE_SIZE) + 1;
	}
	
	if (ydim % STARPU_TILE_SIZE == 0) {
		yt = ydim / STARPU_TILE_SIZE;
	} else {
		yt = (ydim / STARPU_TILE_SIZE) + 1;
	}

#if INTEROP_KERNEL == KERNEL_OMPSS
	interop_ompss_init();

	cpu_set_t cpu_mask;
	sched_getaffinity(0, sizeof(cpu_set_t), &cpu_mask);
	
	int cpu_count = CPU_COUNT(&cpu_mask);
	
	assert((cpu_count % MAX_CPUS_TASK) == 0);
#elif INTEROP_KERNEL == KERNEL_OPENCL
	int num_dev = interop_opencl_init(MAX_CPUS_TASK);
	/*
	 * OpenCL initializes StarpuRM, which creates and sets a new context.
	 * Force running the task in the global context.
	 */
	unsigned ctx = 0;
	starpu_sched_ctx_set_context(&ctx);
#endif
	
	partition_mult_data(xdim, xt, ydim, yt, zdim, A, &A_handle, lda, B, &B_handle, ldb, C, &C_handle, ldc);

	int cur_assigned = 0;
	unsigned x, y;
	for (x = 0; x < xt; x++)
	{
		for (y = 0; y < yt; y++)
		{
			struct starpu_task *task = starpu_task_create();
			struct starpu_dgemm_args *args = malloc(sizeof(*args));
			args->alpha = alpha;
			args->beta = beta;

#if INTEROP_KERNEL == KERNEL_OMPSS
			CPU_ZERO(&args->cpu_mask);

			int i;
			for (i = 0; i < MAX_CPUS_TASK; ++i) {
				CPU_SET(cur_assigned, &args->cpu_mask);
				cur_assigned = (cur_assigned + 1) % cpu_count;
			}
#elif INTEROP_KERNEL == KERNEL_OPENCL
			args->ctx_id = cur_assigned;
			cur_assigned = (cur_assigned + 1) % num_dev;
#endif

			task->cl = &cl;
			task->cl_arg = args;
			task->cl_arg_size = sizeof(*args);
			task->cl_arg_free = 1;

			task->handles[0] = starpu_data_get_sub_data(A_handle, 1, y);
			task->handles[1] = starpu_data_get_sub_data(B_handle, 1, x);
			task->handles[2] = starpu_data_get_sub_data(C_handle, 2, x, y);

			task->flops = 2ULL * (xdim / xt) * (ydim / yt) * zdim;

			ret = starpu_task_submit(task);
			STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");
			starpu_data_wont_use(starpu_data_get_sub_data(C_handle, 2, x, y));
		}
	}

	starpu_task_wait_for_all();

	starpu_data_unpartition(C_handle, STARPU_MAIN_RAM);
	starpu_data_unpartition(B_handle, STARPU_MAIN_RAM);
	starpu_data_unpartition(A_handle, STARPU_MAIN_RAM);

	starpu_data_unregister(A_handle);
	starpu_data_unregister(B_handle);
	starpu_data_unregister(C_handle);

#if INTEROP_KERNEL == KERNEL_OMPSS
	interop_ompss_finalize();
#elif INTEROP_KERNEL == KERNEL_OPENCL
	interop_opencl_finalize();
#endif

	return ret;
}

void interop_starpu_finalize()
{
	starpu_shutdown();
}
