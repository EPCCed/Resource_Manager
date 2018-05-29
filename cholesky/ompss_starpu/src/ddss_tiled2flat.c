#include "../include/ddss.h"

/**
 *
 * 	@file ddss_tiled2flat.c
 *
 * 	@brief LASs-DDSs ddss_tiled2flat routines.
 *
 * 	LASs-DDSs is a software package provided by:
 * 	Barcelona Supercomputing Center - Centro Nacional de Supercomputacion
 *
 * 	@author Pedro Valero-Lara pedro.valero@bsc.es
 * 	@date 2017-05-10
 * 	@reviewer 
 * 	@modified 
 *
 **/

/**
 *  
 *	@ingroup DDSs
 *   
 *  ddss_dsymtiled2flat:
 *	Performs the change of the data layout from tile layout to flat layout 
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
 *
 **/ 

/**
 *
 *	@sa ddss_dscatter_tile
 *	@sa ddss_tile_size
 *
 **/ 

void ddss_dsymtiled2flat( int M, int N, double *A, int LDA, int MT, int NT, 
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
					ddss_dscatter_tile( M, N, 
						&A[m * TILE_SIZE * N + n * TILE_SIZE], LDA, 
						TILE_A[m][n], m, n );
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
					ddss_dscatter_tile( M, N, 
						&A[m * TILE_SIZE * N + n * TILE_SIZE], LDA, 
						TILE_A[m][n], m, n );
				}
			}
		}
	}		

	#pragma oss taskwait	
}
/**
 *  
 *	@ingroup DDSs
 *   
 *  ddss_dscatter_tile:
 *	Performs the copy of a tile from the tile matrix TILE_A to the flat matrix A
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
 *	@sa ddss_dtiled2flat
 *
 **/ 

void ddss_dscatter_tile( int M, int N, double *A, int LDA, 
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
	 		A[i * LDA + j] = TILE_A[i * tile_size_n + j];
	 	}
	}

}
