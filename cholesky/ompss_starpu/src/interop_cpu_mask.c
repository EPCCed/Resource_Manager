#define _GNU_SOURCE
#include <sched.h>
#include <hwloc.h>
#include <hwloc/glibc-sched.h>
#include <stdlib.h>
#include <assert.h>

hwloc_topology_t topology;
cpu_set_t *cpu_masks;

void gen_cpu_masks(int num_kernels, int max_cpus_per_kernel)
{	
	hwloc_topology_init(&topology);
	hwloc_topology_load(topology);
	
	int depth = hwloc_get_type_depth(topology, HWLOC_OBJ_PU);
	assert(depth != HWLOC_TYPE_DEPTH_UNKNOWN);
	
	int num_cpus = hwloc_get_nbobjs_by_depth(topology, depth);
	cpu_masks = malloc(num_kernels * sizeof(*cpu_masks));

	if (max_cpus_per_kernel == -1) {
		max_cpus_per_kernel = num_cpus / num_kernels;
	}
	
	cpu_set_t global_cpu_mask;
	CPU_ZERO(&global_cpu_mask);
	
	/* Assuming we are using all CPUs, and they have consecutive IDs. This is
	 * because in Nanos6, the main function is a task. Therefore, its CPU mask
	 * has only one CPU marked.
	 */
	int cur_cpu = 0;
	int i, j;
	for (i = 0; i < num_kernels; ++i) {
		CPU_ZERO(&cpu_masks[i]);
		
		for (j = cur_cpu; j < cur_cpu + max_cpus_per_kernel; ++j) {
			CPU_SET(j, &cpu_masks[i]);
		}
		
		cur_cpu = cur_cpu + max_cpus_per_kernel;
		assert(cur_cpu <= num_cpus);
	}
}

cpu_set_t *get_cpu_mask(int kernel)
{
	return &cpu_masks[kernel];
}

hwloc_cpuset_t get_hwloc_cpu_mask(int kernel)
{
	hwloc_cpuset_t hwloc_cpuset = hwloc_bitmap_alloc();
	int status = hwloc_cpuset_from_glibc_sched_affinity(topology, hwloc_cpuset, &cpu_masks[kernel], sizeof(cpu_set_t));
    assert(status == 0);

	return hwloc_cpuset;
}
