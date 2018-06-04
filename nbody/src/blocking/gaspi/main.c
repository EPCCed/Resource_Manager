#include "blocking/gaspi/nbody.h"
#include "blocking/gaspi/macros.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <mpi.h>
#include <GASPI.h>
#include <GASPI_Ext.h>

int main(int argc, char** argv)
{
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	assert(provided == MPI_THREAD_MULTIPLE);
	
	gaspi_rank_t rank, rank_size;
#ifdef INTEROPERABILITY
	SUCCESS_OR_DIE(gaspi_proc_init_mode(GASPI_MODE_TASK, GASPI_BLOCK));
#else
	SUCCESS_OR_DIE(gaspi_proc_init(GASPI_BLOCK));
#endif
	SUCCESS_OR_DIE(gaspi_proc_rank(&rank));
	SUCCESS_OR_DIE(gaspi_proc_num(&rank_size));
	
	nbody_conf_t conf = nbody_get_conf(argc, argv);
	assert(conf.num_particles > 0);
	assert(conf.timesteps > 0);
	
	int total_particles = ROUNDUP(conf.num_particles, MIN_PARTICLES);
	int my_particles = total_particles / rank_size;
	assert(my_particles >= BLOCK_SIZE);
	conf.num_particles = my_particles;
	
	conf.num_blocks = my_particles / BLOCK_SIZE;
	assert(conf.num_blocks > 0);
	
	nbody_t nbody = nbody_setup(&conf);
	SUCCESS_OR_DIE(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));
	
	double start = get_time();
	nbody_solve(&nbody, conf.num_blocks, conf.timesteps, conf.time_interval);
	double end = get_time();
	
	nbody_stats(&nbody, &conf, end - start);
	
	if (conf.save_result) nbody_save_particles(&nbody);
	if (conf.check_result) nbody_check(&nbody);
	nbody_free(&nbody);
	
	SUCCESS_OR_DIE(gaspi_proc_term(GASPI_BLOCK));
	
	MPI_Finalize();
	return 0;
}
