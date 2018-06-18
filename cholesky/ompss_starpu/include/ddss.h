/**
 *
 * 	@file ddss.h
 *
 * 	@brief LASs-DDSs header definition.
 *
 * 	BLASs-DDSs is a software package provided by:
 * 	Barcelona Supercomputing Center - Centro Nacional de Supercomputacion
 *
 * 	@author Pedro Valero-Lara pedro.valero@bsc.es
 * 	@date 2017-01-04
 * 	@reviewer 
 * 	@modified 
 *
 **/

#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <assert.h> 
#include <float.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <mkl.h>
#include <mkl_lapacke.h>
#include "ddss_macros.h" 
#include "ddss_types.h"

#ifndef DDSS_H
#define DDSS_H

// LASs-DDSs routines

// LAPACK

int ddss_dpotrf( enum DDSS_UPLO UPLO, int N, double *A, int LDA );

// BLAS

enum DDSS_RETURN kdpotrf( enum DDSS_UPLO UPLO, int N, double *A, int LDA );

// ddss_tile.c routines

int ddss_tile_size( int M, int MT );

// ddss_flat2tiled.c routines

void ddss_dsymtiled2flat( int M, int N, double *A, int LDA, int MT, int NT, 
		double (*TILE_A)[NT][TILE_SIZE * TILE_SIZE], enum DDSS_UPLO UPLO );

void ddss_dgather_tile( int M, int N, double *A, int LDA, 
		double  *TILE_A, int MID, int NID );

// ddss_tiled2flat.c routines

void ddss_dsymflat2tiled( int M, int N, double *A, int LDA, int MT, int NT, 
		double (*TILE_A)[NT][TILE_SIZE * TILE_SIZE], enum DDSS_UPLO UPLO );

void ddss_dscatter_tile( int M, int N, double *A, int LDA, 
		double *TILE_A, int MID, int NID );

#endif
