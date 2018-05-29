#define _GNU_SOURCE
#include <sched.h>
#include <pthread.h>
#include <assert.h>

#include <CL/cl.h>

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
	int *res;
};

cl_device_id *subdevice_id;
cl_context *subdevice_ctx;
int num_dev;

static void offload_opencl_ddss_dgemm(void *void_args)
{
	struct args_dgemm *args = (struct args_dgemm *)void_args;

	ddss_dgemm(args->TRANS_A, args->TRANS_B, args->M, args->N, args->K,
				args->ALPHA, args->A, args->LDA, args->B, args->LDB,
				args->BETA, args->C, args->LDC);
}

int interop_opencl_init(int max_cpus_task)
{
	int ret;
	
	cl_platform_id platform_id;
	clGetPlatformIDs(1, &platform_id, NULL);
	
	cl_device_id device;
	ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_ALL, 1, &device, NULL);
	assert(ret == CL_SUCCESS);

	int cunits;
	ret = clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS,
			sizeof(int), &cunits, NULL);
	assert(ret == CL_SUCCESS);
	
	assert((cunits % max_cpus_task) == 0);
	num_dev = cunits / max_cpus_task;
	
	cl_device_partition_property props[3];
	props[0] = CL_DEVICE_PARTITION_EQUALLY;
	props[1] = max_cpus_task;
	props[2] = 0;
	
	subdevice_id = malloc(num_dev * sizeof(cl_device_id));
	subdevice_ctx = malloc(num_dev * sizeof(cl_context));

	ret = clCreateSubDevices(device, props, num_dev, subdevice_id, NULL);
	assert(ret == CL_SUCCESS);

	int i;
	for (i = 0; i < num_dev; ++i) {
		subdevice_ctx[i] = clCreateContext(0, 1, &subdevice_id[i], NULL, NULL, NULL);
		assert(subdevice_ctx[i] != NULL);
	}
	
	return num_dev;
}

int interop_opencl_dgemm(int ctx_id,
		int TRANS_A, int TRANS_B,
	    int M, int N, int K,
        double ALPHA, double *A, int LDA,
                      double *B, int LDB,
        double BETA,  double *C, int LDC )
{
	assert(ctx_id < num_dev);
	
	cl_command_queue queue = clCreateCommandQueueWithProperties(subdevice_ctx[ctx_id],
																subdevice_id[ctx_id], 0, NULL);
	assert(queue != NULL);

	// OpenCL makes a copy of the given arguments!
	int mem_res;
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
	args.res = &mem_res;
	
	clEnqueueNativeKernelWithType(queue, &offload_opencl_ddss_dgemm, &args, sizeof(args), 0, NULL, NULL,
			CL_NATIVE_KERNEL_OMPSS, 0, NULL, NULL);
	clFinish(queue);
	
	clReleaseCommandQueue(queue);

	return mem_res;
}

void interop_opencl_finalize()
{
	int i;
	for (i = 0; i < num_dev; ++i) {
		clReleaseContext(subdevice_ctx[i]);
	}
}
