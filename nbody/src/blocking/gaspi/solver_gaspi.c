#include "blocking/gaspi/nbody.h"
#include "blocking/gaspi/macros.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>
#include <GASPI.h>
#include <GASPI_Ext.h>

static void calculate_forces(forces_block_t *forces, const particles_block_t *block1, const particles_block_t *block2, const int num_blocks);
static void update_particles(particles_block_t *particles, forces_block_t *forces, const int num_blocks, const float time_interval);
static void exchange_particles(const particles_block_t *sendbuf, particles_block_t *recvbuf, const int num_blocks, const int rank, const int rank_size, nbody_t *nbody);
static void calculate_forces_block(forces_block_t *forces, const particles_block_t *block1, const particles_block_t *block2);
static void update_particles_block(particles_block_t *particles, forces_block_t *forces, const float time_interval);

void nbody_solve(nbody_t *nbody, const int num_blocks, const int timesteps, const float time_interval)
{
	assert(nbody != NULL);
	assert(timesteps > 0);
	
	gaspi_rank_t rank, rank_size;
	SUCCESS_OR_DIE(gaspi_proc_rank(&rank));
	SUCCESS_OR_DIE(gaspi_proc_num(&rank_size));
	
	particles_block_t *local = nbody->local;
	particles_block_t *remote1 = nbody->remote1;
	particles_block_t *remote2 = nbody->remote2;
	forces_block_t *forces = nbody->forces;
	
	for (int t = 0; t < timesteps; t++) {
		particles_block_t *sendbuf = local;
		particles_block_t *recvbuf = remote1;
		
		for (int r = 0; r < rank_size; r++) {
			calculate_forces(forces, local, sendbuf, num_blocks);
			if (r < rank_size - 1) {
				exchange_particles(sendbuf, recvbuf, num_blocks, rank, rank_size, nbody);
			}
			
			particles_block_t *aux = recvbuf;
			recvbuf = (r != 0) ? sendbuf : remote2;
			sendbuf = aux;
		}
		
		update_particles(local, forces, num_blocks, time_interval);
	}
	
	SUCCESS_OR_DIE(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));
}

void calculate_forces_N2(forces_block_t *forces, const particles_block_t *block1, const particles_block_t *block2, const int num_blocks)
{
	for (int i = 0; i < num_blocks; i++) {
		for (int j = 0; j < num_blocks; j++) {
			calculate_forces_block(forces+i, block1+i, block2+j);
		}
	}
}

void calculate_forces_NlogN(forces_block_t *forces, const particles_block_t *block1, const particles_block_t *block2, const int num_blocks)
{
	for (int i = 0; i < num_blocks; i++) {
		for (int j = 0; j < LOG2(num_blocks); j++) {
			calculate_forces_block(forces+i, block1+i, block2+j);
		}
	}
}

void calculate_forces_N(forces_block_t *forces, const particles_block_t *block1, const particles_block_t *block2, const int num_blocks)
{
	for (int i = 0; i < num_blocks - 1; i++) {
		calculate_forces_block(forces+i, block1+i, block2+i+1);
	}
}

void update_particles(particles_block_t *particles, forces_block_t *forces, const int num_blocks, const float time_interval)
{
	for (int i = 0; i < num_blocks; i++) {
		update_particles_block(particles+i, forces+i, time_interval);
	}
}

void exchange_particles(const particles_block_t *sendbuf, particles_block_t *recvbuf, const int num_blocks, const int rank, const int rank_size, nbody_t *nbody)
{
	gaspi_rank_t prev = MOD(rank - 1, rank_size);
	gaspi_rank_t next = MOD(rank + 1, rank_size);
	
	gaspi_segment_id_t send_seg_id = SEGMENT_ID(nbody, sendbuf);
	gaspi_segment_id_t recv_seg_id = SEGMENT_ID(nbody, recvbuf);
	assert(send_seg_id != recv_seg_id);
	
	gaspi_queue_id_t queue;
	SUCCESS_OR_DIE(gaspi_queue_group_get_queue(0, &queue));
	
	const gaspi_notification_id_t recv_not_id = 0;
	const gaspi_notification_id_t send_not_id = 1;
	gaspi_notification_id_t notified_id;
	gaspi_notification_t notified_value;
	
	// Notify the prev rank that the data can be received
	SUCCESS_OR_DIE(gaspi_notify(send_seg_id, prev, send_not_id, 1, queue, GASPI_BLOCK));
	
	// Wait the notification from the next rank to start the transmission
	SUCCESS_OR_DIE(gaspi_notify_waitsome(send_seg_id, send_not_id, 1, &notified_id, GASPI_BLOCK));
	assert(notified_id == send_not_id);
	SUCCESS_OR_DIE(gaspi_notify_reset(send_seg_id, send_not_id, &notified_value));
	assert(notified_value == 1);
	
	// Send the blocks to the next rank
	gaspi_size_t size = num_blocks * sizeof(particles_block_t);
	SUCCESS_OR_DIE(gaspi_write_notify(
			send_seg_id, 0,
			next, recv_seg_id, 0,
			size, recv_not_id, 1,
			queue, GASPI_BLOCK));
	
	// Wait the writes and notification to be sent
	SUCCESS_OR_DIE(gaspi_wait(queue, GASPI_BLOCK));
	
	// Wait the data to be received
	SUCCESS_OR_DIE(gaspi_notify_waitsome(recv_seg_id, recv_not_id, 1, &notified_id, GASPI_BLOCK));
	assert(notified_id == recv_not_id);
	SUCCESS_OR_DIE(gaspi_notify_reset(recv_seg_id, recv_not_id, &notified_value));
	assert(notified_value == 1);
}

void calculate_forces_block(forces_block_t *forces, const particles_block_t *block1, const particles_block_t *block2)
{
	float *x = forces->x;
	float *y = forces->y;
	float *z = forces->z;
	
	const int same_block = (block1 == block2);
	const float *pos_x1 = block1->position_x;
	const float *pos_y1 = block1->position_y;
	const float *pos_z1 = block1->position_z;
	const float *mass1  = block1->mass ;
	
	const float *pos_x2 = block2->position_x;
	const float *pos_y2 = block2->position_y;
	const float *pos_z2 = block2->position_z;
	const float *mass2  = block2->mass;
	
	for (int i = 0; i < BLOCK_SIZE; i++) {
		float fx = x[i], fy = y[i], fz = z[i];
		for (int j = 0; j < BLOCK_SIZE; j++) {
			const float diff_x = pos_x2[j] - pos_x1[i];
			const float diff_y = pos_y2[j] - pos_y1[i];
			const float diff_z = pos_z2[j] - pos_z1[i];
			
			const float distance_squared = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
			const float distance = sqrtf(distance_squared);
			
			float force = 0.0f;
			if (!same_block || distance_squared != 0.0f) {
				force = (mass1[i] / (distance_squared * distance)) * (mass2[j] * gravitational_constant);
			}
			fx += force * diff_x;
			fy += force * diff_y;
			fz += force * diff_z;
		}
		x[i] = fx;
		y[i] = fy;
		z[i] = fz;
	}
}

void update_particles_block(particles_block_t *particles, forces_block_t *forces, const float time_interval)
{
	for (int e = 0; e < BLOCK_SIZE; e++){
		const float mass       = particles->mass[e];
		const float velocity_x = particles->velocity_x[e];
		const float velocity_y = particles->velocity_y[e];
		const float velocity_z = particles->velocity_z[e];
		const float position_x = particles->position_x[e];
		const float position_y = particles->position_y[e];
		const float position_z = particles->position_z[e];
		
		const float time_by_mass       = time_interval / mass;
		const float half_time_interval = 0.5f * time_interval;
		
		const float velocity_change_x = forces->x[e] * time_by_mass;
		const float velocity_change_y = forces->y[e] * time_by_mass;
		const float velocity_change_z = forces->z[e] * time_by_mass;
		const float position_change_x = velocity_x + velocity_change_x * half_time_interval;
		const float position_change_y = velocity_y + velocity_change_y * half_time_interval;
		const float position_change_z = velocity_z + velocity_change_z * half_time_interval;
		
		particles->velocity_x[e] = velocity_x + velocity_change_x;
		particles->velocity_y[e] = velocity_y + velocity_change_y;
		particles->velocity_z[e] = velocity_z + velocity_change_z;
		particles->position_x[e] = position_x + position_change_x;
		particles->position_y[e] = position_y + position_change_y;
		particles->position_z[e] = position_z + position_change_z;
	}
	
	memset(forces, 0, sizeof(forces_block_t));
}

void nbody_stats(const nbody_t *nbody, const nbody_conf_t *conf, double time)
{
	gaspi_rank_t rank, rank_size;
	SUCCESS_OR_DIE(gaspi_proc_rank(&rank));
	SUCCESS_OR_DIE(gaspi_proc_num(&rank_size));
	
	if (!rank) {
		int particles = nbody->num_blocks * BLOCK_SIZE;
		int total_particles = particles * rank_size;
		
		printf("bigo, %s, processes, %d, threads, %d, timesteps, %d, total_particles, %d, particles_per_proc, %d, block_size, %d, blocks_per_proc, %d, time, %.2f, performance, %.2f\n",
			TOSTRING(BIGO), rank_size, 1, nbody->timesteps, total_particles, particles, BLOCK_SIZE,
			nbody->num_blocks, time, nbody_compute_throughput(total_particles, nbody->timesteps, time)
		);
	}
}

