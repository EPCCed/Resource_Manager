#include "blocking/gaspi/nbody.h"
#include "blocking/gaspi/macros.h"

#include <assert.h>
#include <ctype.h>
#include <fcntl.h>
#include <getopt.h>
#include <ieee754.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include <mpi.h>
#include <GASPI.h>
#include <GASPI_Ext.h>

void nbody_generate_particles(const nbody_conf_t *conf, const nbody_file_t *file)
{
	gaspi_rank_t rank_size;
	SUCCESS_OR_DIE(gaspi_proc_num(&rank_size));
	
	char fname[1024];
	sprintf(fname, "%s.in", file->name);
	
	if (!access(fname, F_OK)) {
		return;
	}
	
	struct stat st = {0};
	if (stat("data", &st) == -1) {
		mkdir("data", 0755);
	}
	
	const int fd = open(fname, O_RDWR|O_CREAT|O_TRUNC, S_IRUSR|S_IRGRP|S_IROTH);
	assert(fd >= 0);
	
	const int total_size = file->total_size;
	assert(total_size % PAGE_SIZE == 0);
	
	int err = ftruncate(fd, total_size);
	assert(!err);
	
	particles_block_t * const particles = mmap(NULL, total_size, PROT_WRITE|PROT_READ, MAP_SHARED, fd, 0);
	
	for(int i = 0; i < conf->num_blocks * rank_size; i++) {
		nbody_particle_init(conf, particles+i);
	}
	
	err = munmap(particles, total_size);
	assert(!err);
	
	err = close(fd);
	assert(!err);
}

void nbody_check(const nbody_t *nbody)
{
	gaspi_rank_t rank;
	SUCCESS_OR_DIE(gaspi_proc_rank(&rank));
	
	char fname[1024];
	sprintf(fname, "%s.ref", nbody->file.name);
	if (access(fname, F_OK) != 0) {
		if (!rank) fprintf(stderr, "Warning: %s file does not exist. Skipping the check...\n", fname);
		return;
	}
	
	const int fd = open(fname, O_RDONLY, 0);
	assert(fd >= 0);
	
	particles_block_t *reference = mmap(NULL, nbody->file.size, PROT_READ, MAP_SHARED, fd, nbody->file.offset);
	assert(reference != MAP_FAILED);
	
	int correct, correct_chunk;
	correct_chunk = nbody_compare_particles(nbody->local, reference, nbody->num_blocks);
	
	SUCCESS_OR_DIE(gaspi_allreduce(
			&correct_chunk, &correct, 1,
			GASPI_OP_MIN, GASPI_TYPE_INT,
			GASPI_GROUP_ALL, GASPI_BLOCK));
	
	if (!rank) {
		if (correct) {
			printf("Result validation: OK\n");
		} else {
			printf("Result validation: ERROR\n");
		}
	}
	
	int err = munmap(reference, nbody->file.size);
	assert(!err);
	
	err = close(fd);
	assert(!err);
}

nbody_file_t nbody_setup_file(const nbody_conf_t *conf)
{
	gaspi_rank_t rank, rank_size;
	SUCCESS_OR_DIE(gaspi_proc_rank(&rank));
	SUCCESS_OR_DIE(gaspi_proc_num(&rank_size));
	
	nbody_file_t file;
	file.size = conf->num_blocks * sizeof(particles_block_t);
	file.total_size = rank_size * file.size;
	file.offset = rank * file.size;
	
	sprintf(file.name, "%s-%s-%d-%d-%d", conf->name, TOSTRING(BIGO), rank_size * conf->num_blocks * BLOCK_SIZE, BLOCK_SIZE, conf->timesteps);
	return file;
}

particles_block_t * nbody_load_particles(const nbody_conf_t *conf, const nbody_file_t *file)
{
	char fname[1024];
	sprintf(fname, "%s.in", file->name);
	
	const int fd = open(fname, O_RDONLY, 0);
	assert(fd >= 0);
	
	void * const ptr = mmap(NULL, file->size, PROT_READ|PROT_WRITE, MAP_PRIVATE, fd, file->offset);
	assert(ptr != MAP_FAILED);
	
	int err = close(fd);
	assert(!err);
	
	return ptr;
}

void nbody_setup_gaspi(nbody_t *nbody)
{
	assert(nbody != NULL);
	
	gaspi_info_t *info = &nbody->gaspi_info;
	int particles_size = nbody->num_blocks * sizeof(particles_block_t);
	
	// Commit the default group
	SUCCESS_OR_DIE(gaspi_group_commit(GASPI_GROUP_ALL, GASPI_BLOCK));
	
	// Register a segment for the local buffer of particles
	SUCCESS_OR_DIE(gaspi_segment_use(
			LOCAL_SEGMENT_ID, nbody->local,
			particles_size, GASPI_GROUP_ALL,
			GASPI_BLOCK, 0));
	
	// Register a segment for the 1st remote buffer of particles
	SUCCESS_OR_DIE(gaspi_segment_use(
			REMOTE1_SEGMENT_ID, nbody->remote1,
			particles_size, GASPI_GROUP_ALL,
			GASPI_BLOCK, 0));
	
	// Register a segment for the 2nd remote buffer of particles
	SUCCESS_OR_DIE(gaspi_segment_use(
			REMOTE2_SEGMENT_ID, nbody->remote2,
			particles_size, GASPI_GROUP_ALL,
			GASPI_BLOCK, 0));
	
	// Create the desired number of queues
	SUCCESS_OR_DIE(gaspi_queue_num(&info->precreated_queues));
	SUCCESS_OR_DIE(gaspi_queue_max(&info->max_queues));
	info->max_queues = MIN(info->max_queues, NUM_GASPI_QUEUES);
	assert(info->max_queues > 0);
	
	info->queues = (gaspi_queue_id_t *) malloc(info->max_queues * sizeof(gaspi_queue_id_t));
	assert(info->queues != NULL);
	
	// Some queues are precreated in the GASPI initialization
	for (int i = 0; i < info->precreated_queues; ++i) {
		info->queues[i] = i;
	}
	
	// Create the rest of queues
	for (int i = info->precreated_queues; i < info->max_queues; ++i) {
		SUCCESS_OR_DIE(gaspi_queue_create(&info->queues[i], GASPI_BLOCK));
	}
	
	const gaspi_queue_group_policy_t policy =
#ifdef INTEROPERABILITY
		GASPI_QUEUE_GROUP_POLICY_CPU_RR;
#else
		GASPI_QUEUE_GROUP_POLICY_DEFAULT;
#endif
	
	SUCCESS_OR_DIE(gaspi_queue_group_create(0, 0, info->max_queues, policy));
}

nbody_t nbody_setup(const nbody_conf_t *conf)
{
	gaspi_rank_t rank;
	SUCCESS_OR_DIE(gaspi_proc_rank(&rank));
	
	nbody_file_t file = nbody_setup_file(conf);
	
	if (!rank) nbody_generate_particles(conf, &file);
	SUCCESS_OR_DIE(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));
	
	nbody_t nbody;
	nbody.local = nbody_load_particles(conf, &file);
	assert(nbody.local != NULL);
	
	int particles_size = conf->num_blocks * sizeof(particles_block_t);
	nbody.remote1 = nbody_alloc(particles_size);
	nbody.remote2 = nbody_alloc(particles_size);
	assert(nbody.remote1 != NULL);
	assert(nbody.remote2 != NULL);
	
	memcpy(nbody.remote1, nbody.local, particles_size);
	memcpy(nbody.remote2, nbody.local, particles_size);
	
	nbody.forces = nbody_alloc(conf->num_blocks * sizeof(forces_block_t));
	assert(nbody.forces != NULL);
	
	nbody.num_blocks = conf->num_blocks;
	nbody.timesteps = conf->timesteps;
	nbody.file = file;
	
	nbody_setup_gaspi(&nbody);
	
	return nbody;
}

void nbody_save_particles(const nbody_t *nbody)
{
	SUCCESS_OR_DIE(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));
	
	char fname[1024];
	sprintf(fname, "%s.out", nbody->file.name);
	const int mode = MPI_MODE_CREATE | MPI_MODE_WRONLY;
	
	MPI_File outfile;
	int err = MPI_File_open(MPI_COMM_WORLD, fname, mode, MPI_INFO_NULL, &outfile);
	assert(err == MPI_SUCCESS);
	
	MPI_File_set_view(outfile, nbody->file.offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
	assert(err == MPI_SUCCESS);
	
	MPI_File_write(outfile, nbody->local, nbody->file.size, MPI_BYTE, MPI_STATUS_IGNORE);
	assert(err == MPI_SUCCESS);
	
	MPI_File_close(&outfile);
	assert(err == MPI_SUCCESS);
	
	MPI_Barrier(MPI_COMM_WORLD);
}

void nbody_free_gaspi(nbody_t *nbody)
{
	gaspi_info_t *info = &nbody->gaspi_info;
	
	for (int i = 0; i < info->max_queues; ++i) {
		SUCCESS_OR_DIE(gaspi_wait(info->queues[i], GASPI_BLOCK));
	}
	
	SUCCESS_OR_DIE(gaspi_queue_group_delete(0));
	
	for (int i = info->precreated_queues; i < info->max_queues; ++i) {
		SUCCESS_OR_DIE(gaspi_queue_delete(info->queues[i]));
	}
	
	free(info->queues);
}

void nbody_free(nbody_t *nbody)
{
	const int particles_size = nbody->num_blocks * sizeof(particles_block_t);
	const int forces_size = nbody->num_blocks * sizeof(forces_block_t);
	
	nbody_free_gaspi(nbody);
	
	int err = munmap(nbody->local, particles_size);
	err |= munmap(nbody->remote1, particles_size);
	err |= munmap(nbody->remote2, particles_size);
	err |= munmap(nbody->forces, forces_size);
	assert(!err);
}

