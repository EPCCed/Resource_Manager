#include "../include/ddss.h"

/**
 *
 * 	@file ddss_dpotrf.c
 *
 * 	@brief LASs-DDSs ddss_dpotrf routine.
 *
 * 	LASs-DDSs is a software package provided by:
 * 	Barcelona Supercomputing Center - Centro Nacional de Supercomputacion
 *
 * 	@author Pedro Valero-Lara pedro.valero@bsc.es
 * 	@date 2017-15-08
 * 	@reviewer 
 * 	@modified 
 *
 **/

/**
 *  
 *	@ingroup DDSS
 *   
 *	Performs the Cholesky factorization of a simmetric positive definite 
 *	matrix A:
 *
 *		A = L \times L^T
 *		or
 *		A = U^T \times U
 *
 *	where L is a lower triangular matrix and U is an upper triangular matrix.
 *
**/

/**
 *	
 * 	@param[in]
 *	UPLO	enum DDSS_UPLO.
 *	        UPLO specifies the form of A is stored:
 *       	- Lower: Lower triangle of A is stored. The upper traingular part is
 *       	not referenced.
 *      	- Upper: Upper triangle of A is stored. The lower triangular part is
 *      	not referenced.
 *	
 *	@param[in]
 *	N       int.
 *          N specifies the order of the square matrix A. N >= 0.
 *    
 *	@param[in,out]
 *  A  		double *.
 *	    	A is a pointer to a positive definite matrix of dimension N by LDA.
 *          On exit, if return value is Success, the matrix A is overwriten by 
 *          the factor U or L.
 *   
 *	@param[in]
 *  LDA     int.
 * 	    	LDA specifies the number of columns of A ( row-major order ).
 * 	    	LDA must be at least max( 1, N ).
 *
 **/

/**
 *
 *	@retval Success successful exit
 *	@retval NoSuccess unsuccessful exit
 *
 **/

/**
 *
 *	@sa kdpotrf
 *
 **/ 

int ddss_dpotrf( enum DDSS_UPLO UPLO, int N, double *A, int LDA )
{

	// Argument checking
	if ( ( UPLO != Upper ) && ( UPLO != Lower ) )
	{
		fprintf( stderr, "Illegal value of UPLO, in ddss_dportf code\n" );
		return NoSuccess;
	}
	
	if ( N < 0 )
	{
		fprintf( stderr, "Illegal value of N, in ddss_dportf code\n" );
		return NoSuccess;
	}

	if ( LDA < MAX( 1, N ) )
	{
		fprintf( stderr, "Illegal value of LDA, in ddss_dportf code\n" );
		return NoSuccess;
	}

	// Quick return
	if ( MAX( N, 0 ) == 0 )
	{
		return Success;
	}
	
	return kdpotrf( UPLO, N, A, LDA );
	
}


