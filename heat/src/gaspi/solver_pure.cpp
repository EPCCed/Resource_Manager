#include <mpi.h>
#include <GASPI.h>
#include <GASPI_Ext.h>
#include <algorithm>
#include <cassert>

#include "common/heat.hpp"
#include "gaspi/macros.hpp"
#include "gaspi/utils.hpp"


inline void solveBlock(block_t *matrix, int nbx, int nby, int bx, int by)
{
	block_t &targetBlock = matrix[bx*nby + by];
	const block_t &centerBlock = matrix[bx*nby + by];
	const block_t &topBlock    = matrix[(bx-1)*nby + by];
	const block_t &leftBlock   = matrix[bx*nby + (by-1)];
	const block_t &rightBlock  = matrix[bx*nby + (by+1)];
	const block_t &bottomBlock = matrix[(bx+1)*nby + by];
	
	double sum = 0.0;
	for (int x = 0; x < BSX; ++x) {
		const row_t &topRow = (x > 0) ? centerBlock[x-1] : topBlock[BSX-1];
		const row_t &bottomRow = (x < BSX-1) ? centerBlock[x+1] : bottomBlock[0];
		
		for (int y = 0; y < BSY; ++y) {
			double left = (y > 0) ? centerBlock[x][y-1] : leftBlock[x][BSY-1];
			double right = (y < BSY-1) ? centerBlock[x][y+1] : rightBlock[x][0];
			
			double value = 0.25 * (topRow[y] + bottomRow[y] + left + right);
			double diff = value - targetBlock[x][y];
			sum += diff * diff;
			targetBlock[x][y] = value;
		}
	}
}

inline void sendFirstComputeRow(block_t *matrix, int nbx, int nby, int rank, int rank_size, gaspi_info_t &info)
{
	const gaspi_offset_t baseOffset = FCR_OFFSET(nbx, nby);
	const gaspi_offset_t remoteBaseOffset = LB_OFFSET(nbx, nby);
	
	gaspi_queue_id_t queue;
	SUCCESS_OR_DIE(gaspi_queue_group_get_queue(0, &queue));
	
	WAIT_FOR_ENTRIES(queue, (nby-2)*2);
	
	for (int by = 1; by < nby-1; ++by) {
		SUCCESS_OR_DIE(gaspi_write_notify(
				SEGMENT_ID, OFFSET(baseOffset, by),
				rank - 1, SEGMENT_ID, OFFSET(remoteBaseOffset, by),
				BSY * sizeof(double), nby + by, 1,
				queue, GASPI_BLOCK));
	}
}

inline void sendLastComputeRow(block_t *matrix, int nbx, int nby, int rank, int rank_size, gaspi_info_t &info)
{
	const gaspi_offset_t baseOffset = LCR_OFFSET(nbx, nby);
	const gaspi_offset_t remoteBaseOffset = UB_OFFSET(nbx, nby);
	
	gaspi_queue_id_t queue;
	SUCCESS_OR_DIE(gaspi_queue_group_get_queue(0, &queue));
	
	for (int by = 1; by < nby-1; ++by) {
		SUCCESS_OR_DIE(gaspi_write_notify(
				SEGMENT_ID, OFFSET(baseOffset, by),
				rank + 1, SEGMENT_ID, OFFSET(remoteBaseOffset, by),
				BSY * sizeof(double), by, 1,
				queue, GASPI_BLOCK));
	}
	
	SUCCESS_OR_DIE(gaspi_wait(queue, GASPI_BLOCK));
}

inline void receiveUpperAndLowerBorders(block_t *matrix, int nbx, int nby, int rank, int rank_size, gaspi_info_t &info)
{
	gaspi_notification_id_t notifiedId;
	gaspi_notification_t notifiedValue;
	
	int expectedNotifications = (rank > 0 && rank < rank_size - 1) ? 2 * (nby - 2) : nby - 2;
	
	while (expectedNotifications) {
		SUCCESS_OR_DIE(gaspi_notify_waitsome(
				SEGMENT_ID, 0, 2 * nby,
				&notifiedId, GASPI_BLOCK));
		
		SUCCESS_OR_DIE(gaspi_notify_reset(
				SEGMENT_ID, notifiedId,
				&notifiedValue));
		assert(notifiedValue == 1);
		
		--expectedNotifications;
	}
}

inline void solveGaussSeidel(block_t *matrix, int nbx, int nby, int rank, int rank_size, gaspi_info_t &info)
{
	if (rank != 0) {
		sendFirstComputeRow(matrix, nbx, nby, rank, rank_size, info);
	}
	
	if (rank_size > 1) {
		receiveUpperAndLowerBorders(matrix, nbx, nby, rank, rank_size, info);
	}
	
	for (int bx = 1; bx < nbx-1; ++bx) {
		for (int by = 1; by < nby-1; ++by) {
			solveBlock(matrix, nbx, nby, bx, by);
		}
	}
	
	if (rank != rank_size - 1) {
		sendLastComputeRow(matrix, nbx, nby, rank, rank_size, info);
	}
}

template<>
double solve(block_t *matrix, int rowBlocks, int colBlocks, int timesteps, gaspi_info_t &info)
{
	gaspi_rank_t rank, rank_size;
	SUCCESS_OR_DIE(gaspi_proc_rank(&rank));
	SUCCESS_OR_DIE(gaspi_proc_num(&rank_size));
	
	for (int t = 0; t < timesteps; ++t) {
		solveGaussSeidel(matrix, rowBlocks, colBlocks, rank, rank_size, info);
	}
	
	SUCCESS_OR_DIE(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));
	
	return 0.0;
}

