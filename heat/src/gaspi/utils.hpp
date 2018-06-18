#ifndef UTILS_HPP
#define UTILS_HPP

#include <GASPI.h>

#include "common/heat.hpp"

#define SEGMENT_ID 0

#ifdef INTEROPERABILITY
#define NUM_GASPI_QUEUES 64
#else
#define NUM_GASPI_QUEUES 1
#endif

struct gaspi_info_t {
	gaspi_queue_id_t *queues;
	gaspi_number_t precreatedQueues;
	gaspi_number_t maxQueues;
	
	gaspi_info_t() :
		queues(nullptr),
		precreatedQueues(0),
		maxQueues(0)
	{
	}
};

gaspi_info_t setupGaspiInfo(const HeatConfiguration &conf, int blocksPerRank);
void freeGaspiInfo(gaspi_info_t &info);

#endif // UTILS_HPP

