#ifndef NBODY_GASPI_H
#define NBODY_GASPI_H

#include "blocking/common/nbody.h"

#include <GASPI.h>

#ifdef INTEROPERABILITY
#define NUM_GASPI_QUEUES 64
#else
#define NUM_GASPI_QUEUES 1
#endif

// Application structures
struct nbody_file_t {
	size_t total_size;
	size_t size;
	size_t offset;
	char name[1000];
};

typedef struct {
	gaspi_queue_id_t *queues;
	gaspi_number_t precreated_queues;
	gaspi_number_t max_queues;
	gaspi_number_t max_requests;
	gaspi_number_t inflight_requests;
} gaspi_info_t;

struct nbody_t {
	particles_block_t *local;
	particles_block_t *remote1;
	particles_block_t *remote2;
	forces_block_t *forces;
	int num_blocks;
	int timesteps;
	nbody_file_t file;
	gaspi_info_t gaspi_info;
};

#endif // NBODY_GASPI_H

