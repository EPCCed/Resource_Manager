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
#if defined(LASs_WITH_MKL)
#include <mkl.h>
#include <mkl_lapacke.h>
#else
#include <cblas.h>
#include <lapacke.h>
#endif
#include "ddss_macros.h" 
#include "ddss_types.h"

#ifndef DDSS_H
#define DDSS_H

// LASs-DDSs routines

// BLAS-3

int ddss_dgemm( enum DDSS_TRANS TRANS_A, enum DDSS_TRANS TRANS_B,
	    int M, int N, int K,
        double ALPHA, double *A, int LDA,
                      double *B, int LDB,
        double BETA,  double *C, int LDC );

int ddss_dsymm( enum DDSS_SIDE SIDE, enum DDSS_UPLO UPLO, 
		int M, int N,
		double ALPHA, double *A, int LDA,
					  double *B, int LDB,
		double BETA,  double *C, int LDC ); 

int ddss_dtrsm( enum DDSS_SIDE SIDE, enum DDSS_UPLO UPLO, 
		enum DDSS_TRANS TRANS_A,
        enum DDSS_DIAG DIAG, int M, int N,
        const double ALPHA, double *A, int LDA,
                            double *B, int LDB);

/*
int ddss_dsyrk( const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, 
		MKL_INT N, MKL_INT K,
        const double ALPHA, double *A, MKL_INT LDA,
    	const double BETA,  double *C, MKL_INT LDC);
*/

// LAPACK

int ddss_dnpgetrf( int M, int N, double *A, int LDA );

int ddss_dpotrf( enum DDSS_UPLO UPLO, int N, double *A, int LDA );

int ddss_dgeqrf( int M, int N, double *A, int LDA, double *T );

// BLAS

enum DDSS_RETURN kdgemm( enum DDSS_TRANS TRANS_A, enum DDSS_TRANS TRANS_B,
		int M, int N, int K,	
		const double ALPHA, double *A, int LDA,
							double *B, int LDB,
		const double BETA,  double *C, int LDC );

enum DDSS_RETURN kdsymm( enum DDSS_SIDE SIDE, enum DDSS_UPLO UPLO, 
		int M, int N,
		const double ALPHA, double *A, int LDA,
					  		double *B, int LDB,
		const double BETA,  double *C, int LDC );

enum DDSS_RETURN kdtrsm( enum DDSS_SIDE SIDE, enum DDSS_UPLO UPLO, 
		enum DDSS_TRANS TRANS_A,
        enum DDSS_DIAG DIAG, int M, int N,
        const double ALPHA, double* A, int LDA,
                            double *B, int LDB);

/*
enum DDSS_RETURN kdsyrk( const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, 
		MKL_INT N, MKL_INT K,
        const double ALPHA, double* A, MKL_INT LDA,
        const double BETA, double *C, MKL_INT LDC);
*/

// LAPACK

enum DDSS_RETURN kdnpgetrf( int M, int N, double *A, int LDA );

enum DDSS_RETURN kdpotrf( enum DDSS_UPLO UPLO, int N, double *A, int LDA );

enum DDSS_RETURN kgeqrf( int M, int N, double *A, int LDA, double *T );

enum DDSS_RETURN dnpgetrf( int M, int N, double *A, int LDA );

// ddss_tile.c routines

int ddss_tile_size( int M, int MT );

// ddss_flat2tiled.c routines

void ddss_dflat2tiled( int M, int N, double *A, int LDA, int MT, int NT, 
		double (*TILE_A)[NT][TILE_SIZE * TILE_SIZE] );

void ddss_dsymtiled2flat( int M, int N, double *A, int LDA, int MT, int NT, 
		double (*TILE_A)[NT][TILE_SIZE * TILE_SIZE], enum DDSS_UPLO UPLO );

void ddss_dgather_tile( int M, int N, double *A, int LDA, 
		double  *TILE_A, int MID, int NID );

// ddss_tiled2flat.c routines

void ddss_dtiled2flat( int M, int N, double *A, int LDA, int MT, int NT, 
		double (*TILE_A)[NT][TILE_SIZE * TILE_SIZE] );

void ddss_dsymflat2tiled( int M, int N, double *A, int LDA, int MT, int NT, 
		double (*TILE_A)[NT][TILE_SIZE * TILE_SIZE], enum DDSS_UPLO UPLO );

void ddss_dscatter_tile( int M, int N, double *A, int LDA, 
		double *TILE_A, int MID, int NID );

#endif
