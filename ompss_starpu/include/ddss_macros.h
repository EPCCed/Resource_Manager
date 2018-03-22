/**
 *
 * 	@file ddss_macros.h
 *
 * 	@brief Macros definition.
 *
 * 	LASs-DDSs is a software package provided by:
 * 	Barcelona Supercomputing Center - Centro Nacional de Supercomputacion
 *
 * 	@author Pedro Valero-Lara pedro.valero@bsc.es
 * 	@date 2017-01-02
 * 	@reviewer 
 * 	@modified 
 *
 **/

#ifndef DDSS_MACROS_H
#define DDSS_MACROS_H

// VERBOSE (activate prints in testings)
//#define VERBOSE

// Tile size
#define TILE_SIZE 4096

// Return the max/min of two numbers
#define MAX( a, b ) ( ( ( a ) > ( b ) ) ? ( a ) : ( b ) )
#define MIN( a, b ) ( ( ( a ) < ( b ) ) ? ( a ) : ( b ) )

// Number of operations 
// POTRF
// Number of multiplications in POTRF
#define FMULS_POTRF(n_) ( (1./6.) * (n_) * (n_) * (n_) + (0.5) * (n_) * (n_) \
	+ (1./3.) * n_ )
// Number of additions in POTRF
#define FADDS_POTRF(n_) ( (1./6.) * (n_) * (n_) * (n_) - (1./6.) * (n_))
// Flops in DPOTRF
#define FLOPS_DPOTRF(n_) ( FMULS_POTRF((double)(n_)) + FADDS_POTRF((double)(n_)))

#endif
