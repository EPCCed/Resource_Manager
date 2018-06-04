#include <GASPI.h>
#include <GASPI_Ext.h>
#include <algorithm>
#include <cassert>

#include "common/heat.hpp"
#include "gaspi/macros.hpp"
#include "gaspi/utils.hpp"

gaspi_info_t setupGaspiInfo(const HeatConfiguration &conf, int blocksPerRank)
{
	gaspi_info_t info;
	
	// Commit the default group
	SUCCESS_OR_DIE(gaspi_group_commit(GASPI_GROUP_ALL, GASPI_BLOCK));
	
	// Register a segment for the local buffer of particles
	SUCCESS_OR_DIE(gaspi_segment_use(
			SEGMENT_ID, conf.matrix,
			blocksPerRank * sizeof(block_t),
			GASPI_GROUP_ALL, GASPI_BLOCK, 0));
	
	// Create the desired number of queues
	SUCCESS_OR_DIE(gaspi_queue_num(&info.precreatedQueues));
	SUCCESS_OR_DIE(gaspi_queue_max(&info.maxQueues));
	info.maxQueues = std::min((int)info.maxQueues, NUM_GASPI_QUEUES);
	assert(info.maxQueues > 0);
	
	info.queues = (gaspi_queue_id_t *) malloc(info.maxQueues * sizeof(gaspi_queue_id_t));
	assert(info.queues != NULL);
	
	// Some queues are precreated in the GASPI initialization
	for (int i = 0; i < info.precreatedQueues; ++i) {
		info.queues[i] = i;
	}
	
	// Create the rest of queues
	for (int i = info.precreatedQueues; i < info.maxQueues; ++i) {
		SUCCESS_OR_DIE(gaspi_queue_create(&info.queues[i], GASPI_BLOCK));
	}
	
#ifdef INTEROPERABILITY
	const gaspi_queue_group_policy_t policy =
		GASPI_QUEUE_GROUP_POLICY_CPU_RR;
#else
	const gaspi_queue_group_policy_t policy =
		GASPI_QUEUE_GROUP_POLICY_DEFAULT;
#endif
	
	// Create the queue group with the previous queues
	SUCCESS_OR_DIE(gaspi_queue_group_create(0, 0, info.maxQueues, policy));
	
	return info;
}

void freeGaspiInfo(gaspi_info_t &info)
{
	for (int i = 0; i < info.maxQueues; ++i) {
		SUCCESS_OR_DIE(gaspi_wait(info.queues[i], GASPI_BLOCK));
	}
	
	SUCCESS_OR_DIE(gaspi_queue_group_delete(0));
	
	for (int i = info.precreatedQueues; i < info.maxQueues; ++i) {
		SUCCESS_OR_DIE(gaspi_queue_delete(info.queues[i]));
	}
	
	free(info.queues);
}

