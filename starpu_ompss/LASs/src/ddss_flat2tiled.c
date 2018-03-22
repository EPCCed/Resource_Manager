#include "../include/ddss.h"

/**
 *
 * 	@file ddss_flat2tiled.c
 *
 * 	@brief LASs-DDSs ddss_flat2tiled routines.
 *
 * 	LASs-DDSs is a software package provided by:
 * 	Barcelona Supercomputing Center - Centro Nacional de Supercomputacion
 *
 * 	@author Pedro Valero-Lara pedro.valero@bsc.es
 * 	@date 2017-05-11
 * 	@reviewer 
 * 	@modified 
 *
 **/

/**
 *  
 *	@ingroup DDSs
 *   
 *  ddss_dflat2tiled:
 *	Performs the change of the data layout from flat layout to tiled layout 
 *  according to row-major order.
 *
**/

/**
 * 
 * 	@param[in]
 *	M		int.	
 *	        M specifies the number of rows of the flat matrix.
 * 
 * 	@param[in]
 *	N		int.	
 *	        N specifies the number of columns of the flat matrix.
 *	
 * 	@param[in]
 *	A		double *.	
 *	        A is a pointer to the flat matrix.
 *
 * 	@param[in]
 *	LDA		int.
 *  		LDA specifies the number of columns ( row-major order ) of matrix A.
 *
 * 	@param[in]
 *	MT		int.	
 *	        MT specifies the number of rows of the matrix TILE_A.
 *			
 * 	@param[in]
 *	NT		int.
 *  		NT specifies the number of columns of the matrix TILE_A.
 *
 * 	@param[in]
 *	TILE_A	double *.
 *  		TILE_A is a pointer to the tile matrix.
 *
 **/ 

/**
 *
 *	@sa ddss_dgather_tile
 *	@sa ddss_tile_size
 *
 **/ 

void ddss_dflat2tiled( int M, int N, double *A, int LDA, int MT, int NT, 
		double (*TILE_A)[NT][TILE_SIZE * TILE_SIZE] )
{
	
	// Local variables
	int m, n;
  	
	for ( m = 0; m < MT; m++ )
	{
		for ( n = 0; n < NT; n++ )
		{
			#pragma oss task inout(TILE_A[m][n])
			{ 
				ddss_dgather_tile( M, N, &A[m * TILE_SIZE * N + n * TILE_SIZE], 
					LDA, TILE_A[m][n], m, n );
			}
		}
	}

}

/**
 *  
 *	@ingroup DDSs
 *   
 *  ddss_dsymflat2tiled:
 *	Performs the change of the data layout from flat layout to tiled layout 
 *  for symmetric matrices according to row-major order.
 *
**/

/**
 * 
 * 	@param[in]
 *	M		int.	
 *	        M specifies the number of rows of the flat matrix.
 * 
 * 	@param[in]
 *	N		int.	
 *	        N specifies the number of columns of the flat matrix.
 *	
 * 	@param[in]
 *	A		double *.	
 *	        A is a pointer to the flat matrix.
 *
 * 	@param[in]
 *	LDA		int.
 *  		LDA specifies the number of columns ( row-major order ) of matrix A.
 *
 * 	@param[in]
 *	MT		int.	
 *	        MT specifies the number of rows of the matrix TILE_A.
 *			
 * 	@param[in]
 *	NT		int.
 *  		NT specifies the number of columns of the matrix TILE_A.
 *
 * 	@param[in]
 *	TILE_A	double *.
 *  		TILE_A is a pointer to the tile matrix.
 *
 * 	@param[in]
 *	UPLO	enum DDSS_UPLO.
 *	        UPLO specifies the form of A is stored:
 *       	- Lower: Lower triangle of A is stored. The upper traingular part is
 *       	not referenced.
 *      	- Upper: Upper triangle of A is stored. The lower triangular part is
 *      	not referenced.

 **/ 

/**
 *
 *	@sa ddss_dgather_tile
 *	@sa ddss_tile_size
 *
 **/ 

void ddss_dsymflat2tiled( int M, int N, double *A, int LDA, int MT, int NT, 
		double (*TILE_A)[NT][TILE_SIZE * TILE_SIZE], enum DDSS_UPLO UPLO )
{
	
	// Local variables
	int m, n;

	if ( UPLO == Upper )
	{
		for ( m = 0; m < MT; m++ )
		{
			for ( n = NT-1; n >= m; n-- )
			{
				#pragma oss task inout(TILE_A[m][n])
				{ 
					ddss_dgather_tile( M, N, 
						&A[m * TILE_SIZE * N + n * TILE_SIZE], 
						LDA, TILE_A[m][n], m, n );
				}
			}
		}
	}
	else if ( UPLO == Lower )
	{
		for ( m = 0; m < MT; m++ )
		{
			for ( n = 0; n <= m; n++ )
			{
				#pragma oss task inout(TILE_A[m][n])
				{ 
					ddss_dgather_tile( M, N, 
						&A[m * TILE_SIZE * N + n * TILE_SIZE], 
						LDA, TILE_A[m][n], m, n );
				}
			}
		}
	}

}


/**
 *  
 *	@ingroup DDSs
 *   
 *  ddss_dgather_tile:
 *	Performs the copy of a tile from the flat matrix A to the tile matrix TILE_A
 * 	for the MT, NT tile.
 *
**/

/**
 *
 * 	@param[in]
 *	M		int.	
 *	        M specifies the number of rows of the flat matrix.
 * 
 * 	@param[in]
 *	N		int.	
 *	        N specifies the number of columns of the flat matrix.
 *	
 * 	@param[in]
 *	A		double *.	
 *	        A is a pointer to the flat matrix.
 *			
 * 	@param[in]
 *	LDA		int.
 *  		LDA specifies the number of columns ( row-major order ) of matrix A.
 *
 * 	@param[in]
 *	TILE_A	double *.
 *  		TILE_A is a pointer to the tile matrix.
 *
 * 	@param[in]
 *	MID		int.
 *  		MID specifies the row id of the tile.
 *
 * 	@param[in]
 *	NID		int.
 *  		NID specifies the column id of the tile. 
 **/

/**
 *
 *	@sa ddss_tile_size
 *	@sa ddss_dflat2tiled
 *
 **/ 

void ddss_dgather_tile( int M, int N, double *A, int LDA, 
		double *TILE_A, int MID, int NID )
{

	//Local variables
	int i, j;
	int tile_size_m, tile_size_n;
	 
	tile_size_m = ddss_tile_size( M, MID );
	tile_size_n = ddss_tile_size( N, NID );

	for ( i = 0; i < tile_size_m; i++ )
	{
		for ( j = 0; j < tile_size_n; j++ )
		{
			TILE_A[i * tile_size_n + j] = A[i * LDA + j];
		}
	}

}
