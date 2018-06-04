#include <mpi.h>
#include <GASPI.h>
#include <GASPI_Ext.h>
#include <algorithm>
#include <cassert>

#include "common/heat.hpp"
#include "gaspi/macros.hpp"
#include "gaspi/utils.hpp"

#ifdef INTEROPERABILITY
#define BLOCK_MODE GASPI_BLOCK_TASK
#else
#define BLOCK_MODE GASPI_BLOCK
#endif

#ifdef INTEROPERABILITY
int *serial = nullptr;
#else
int *serial = (int *) 1;
#endif

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
	
	for (int by = 1; by < nby-1; ++by) {
		#pragma oss task label(send first compute row) in(([nbx][nby]matrix)[1][by]) inout(*serial)
		{
			gaspi_queue_id_t queue;
			SUCCESS_OR_DIE(gaspi_queue_group_get_queue(0, &queue));
			
			SUCCESS_OR_DIE(gaspi_write_notify(
					SEGMENT_ID, OFFSET(baseOffset, by),
					rank - 1, SEGMENT_ID, OFFSET(remoteBaseOffset, by),
					BSY * sizeof(double), nby + by, 1,
					queue, GASPI_BLOCK));
			
			SUCCESS_OR_DIE(gaspi_wait(queue, BLOCK_MODE));
		}
	}
}

inline void sendLastComputeRow(block_t *matrix, int nbx, int nby, int rank, int rank_size, gaspi_info_t &info)
{
	const gaspi_offset_t baseOffset = LCR_OFFSET(nbx, nby);
	const gaspi_offset_t remoteBaseOffset = UB_OFFSET(nbx, nby);
	
	for (int by = 1; by < nby-1; ++by) {
		#pragma oss task label(send last compute row) in(([nbx][nby]matrix)[nbx-2][by]) inout(*serial)
		{
			gaspi_queue_id_t queue;
			SUCCESS_OR_DIE(gaspi_queue_group_get_queue(0, &queue));
			
			SUCCESS_OR_DIE(gaspi_write_notify(
					SEGMENT_ID, OFFSET(baseOffset, by),
					rank + 1, SEGMENT_ID, OFFSET(remoteBaseOffset, by),
					BSY * sizeof(double), by, 1,
					queue, GASPI_BLOCK));
			
			SUCCESS_OR_DIE(gaspi_wait(queue, BLOCK_MODE));
		}
	}
}

inline void receiveUpperBorder(block_t *matrix, int nbx, int nby, int rank, int rank_size, gaspi_info_t &info)
{
	gaspi_notification_id_t notifiedId;
	gaspi_notification_t notifiedValue;
	
	for (int by = 1; by < nby-1; ++by) {
		#pragma oss task label(receive upper border) out(([nbx][nby]matrix)[0][by]) inout(*serial)
		{
			SUCCESS_OR_DIE(gaspi_notify_waitsome(
					SEGMENT_ID, by, 1,
					&notifiedId, BLOCK_MODE));
			assert(notifiedId == by);
			
			SUCCESS_OR_DIE(gaspi_notify_reset(
					SEGMENT_ID, notifiedId,
					&notifiedValue));
			assert(notifiedValue == 1);
		}
	}
}

inline void receiveLowerBorder(block_t *matrix, int nbx, int nby, int rank, int rank_size, gaspi_info_t &info)
{
	gaspi_notification_id_t notifiedId;
	gaspi_notification_t notifiedValue;
	
	for (int by = 1; by < nby-1; ++by) {
		#pragma oss task label(receive lower border) out(([nbx][nby]matrix)[nbx-1][by]) inout(*serial)
		{
			gaspi_notification_id_t notifiedId;
			gaspi_notification_t notifiedValue;
			
			SUCCESS_OR_DIE(gaspi_notify_waitsome(
					SEGMENT_ID, nby + by, 1,
					&notifiedId, BLOCK_MODE));
			assert(notifiedId == nby + by);
		
			SUCCESS_OR_DIE(gaspi_notify_reset(
					SEGMENT_ID, notifiedId,
					&notifiedValue));
			assert(notifiedValue == 1);
		}
	}
}

inline void solveGaussSeidel(block_t *matrix, int nbx, int nby, int rank, int rank_size, gaspi_info_t &info)
{
	if (rank != 0) {
		sendFirstComputeRow(matrix, nbx, nby, rank, rank_size, info);
		receiveUpperBorder(matrix, nbx, nby, rank, rank_size, info);
	}
	
	if (rank != rank_size - 1) {
		receiveLowerBorder(matrix, nbx, nby, rank, rank_size, info);
	}
	
	for (int bx = 1; bx < nbx-1; ++bx) {
		for (int by = 1; by < nby-1; ++by) {
			#pragma oss task label(gauss seidel)     \
					in(([nbx][nby]matrix)[bx-1][by]) \
					in(([nbx][nby]matrix)[bx][by-1]) \
					in(([nbx][nby]matrix)[bx][by+1]) \
					in(([nbx][nby]matrix)[bx+1][by]) \
					inout(([nbx][nby]matrix)[bx][by])
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
	#pragma oss taskwait
	
	SUCCESS_OR_DIE(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));
	
	return 0.0;
}

