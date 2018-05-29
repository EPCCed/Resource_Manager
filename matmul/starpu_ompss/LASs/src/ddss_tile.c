#include "../include/ddss.h"

/**
 *
 * 	@file ddss_tile.c
 *
 * 	@brief LASs-DDSs ddss_tile routines.
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
 *  tile_size:
 *	Computes the size of the tile passed as parameter.
 *
**/

/**
 *	
 * 	@param[in]
 *	M		int.	
 *	        M specifies the size ( rows of columns ) of the matrix.
 *			
 * 	@param[in]
 *	MT		int.
 *  		MT specifies the id of the tile.
 *
 **/

/**
 *
 *	@retval int size of the tile passed as parameter.
 *
 **/

/**
 *
 *	@sa ddss_gather_tile
 *
 **/ 

int ddss_tile_size( int M, int MT )
{

	if ( M - ( MT * TILE_SIZE ) > TILE_SIZE )
		return TILE_SIZE;
	else
		return M - ( MT * TILE_SIZE );

}
