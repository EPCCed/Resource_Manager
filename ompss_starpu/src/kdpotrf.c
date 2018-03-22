#include "../include/ddss.h"
#include "../include/interop.h"

/**
 *
 * 	@file kdpotrf.c
 *
 * 	@brief LASs-DDSs kdgemm routine.
 *
 * 	LASs-DDSs is a software package provided by:
 * 	Barcelona Supercomputing Center - Centro Nacional de Supercomputacion
 *
 * 	@author Pedro Valero-Lara pedro.valero@bsc.es
 * 	@date 2017-24-10
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
 *	    	A is a pointer to a positive definite matrix of dimension LDA by N.
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
 *	@retval Success sucessful exit
 *	@retval NoSuccess unsucessful exit
 *
 **/

/**
 *
 *	@sa ddss_dpotrf
 *	@sa ddss_tile
 *	@sa ddss_symflat2tiled
 *	@sa ddss_symtiled2flat
 *
 **/ 

enum DDSS_RETURN kdpotrf( enum DDSS_UPLO UPLO, int N, double *A, int LDA )
{

	// Local variables
	int nt;
	int mi, mmi, ki, ni;
	int Am;
	int An;
	int tile_size_m;
	int tile_size_mm;
	int tile_size_n;
	int tile_size_k;
	int k_check;
	double betat;
	enum DDSS_RETURN ret;

	// Number of tiles
	if ( N % TILE_SIZE == 0 )
	{
		nt = N / TILE_SIZE;
	}
	else
	{
		nt = ( N / TILE_SIZE ) + 1;
	}
	
    /****************************
	--Tile A matrix declaration--
	****************************/

	Am = nt;
	An = nt;
	
	/***************************
	--Tile A matrix allocation--
	***************************/
		
	double (*TILE_A)[An][TILE_SIZE * TILE_SIZE] = malloc ( Am * An * 
		TILE_SIZE * TILE_SIZE * sizeof( double ) );
	
	if ( TILE_A == NULL)
	{
		fprintf( stderr, "Failure in kdpotrf for matrix TILE_A\n" );
		return NoSuccess;
	}

	/************************************************
	// --From flat data layout to tiled data layout--	
	************************************************/

	ddss_dsymflat2tiled( N, N, A, LDA, nt, nt, TILE_A, UPLO );

	/**************
	--DPOTRF tile--
	**************/	

	if ( UPLO == Lower )
	{ 	
		// --UPLO = Lower--	
		for ( ki = 0; ki < nt; ki++ )
		{
			tile_size_k = ddss_tile_size( N, ki );
	
			#pragma oss task inout( TILE_A[ki][ki] ) \
				shared( TILE_A ) \
				firstprivate( ki )
#ifdef LASs_WITH_CHAMELEON
			interop_dpotrf(UPLO, tile_size_k, TILE_A[ki][ki], 
							 tile_size_k );
#else
			LAPACKE_dpotrf_work( LAPACK_ROW_MAJOR, 
							 'L',
        					 tile_size_k, TILE_A[ki][ki], 
							 tile_size_k );
#endif

			for ( mi = ki + 1; mi < nt; mi++ )
			{
				tile_size_m = ddss_tile_size( N, mi );
			
				#pragma oss task in( TILE_A[ki][ki] ) \
					inout(TILE_A[mi][ki]) \
					shared( TILE_A ) \
					firstprivate( mi, ki )
#ifdef LASs_WITH_CHAMELEON
				interop_dtrsm(( CBLAS_SIDE ) Right, ( CBLAS_UPLO ) Lower, 
					     	 ( CBLAS_TRANSPOSE ) Trans, ( CBLAS_DIAG ) NonUnit, 
						 	 tile_size_m, TILE_SIZE, 
						 	 1.0, TILE_A[ki][ki], TILE_SIZE, 
							  	  TILE_A[mi][ki], TILE_SIZE );
#else
				cblas_dtrsm( CblasRowMajor, 
					     	 ( CBLAS_SIDE ) Right, ( CBLAS_UPLO ) Lower, 
					     	 ( CBLAS_TRANSPOSE ) Trans, ( CBLAS_DIAG ) NonUnit, 
						 	 tile_size_m, TILE_SIZE, 
						 	 1.0, TILE_A[ki][ki], TILE_SIZE, 
							  	  TILE_A[mi][ki], TILE_SIZE );
#endif
			}
			for ( mmi = ki + 1; mmi < nt; mmi++ )
			{
				tile_size_mm = ddss_tile_size( N, mmi );
			
				#pragma oss task in( TILE_A[mmi][ki] ) \
					inout( TILE_A[mmi][mmi] ) \
					shared( TILE_A ) \
					firstprivate( mmi, ki )
#ifdef LASs_WITH_CHAMELEON
				interop_dsyrk(( CBLAS_UPLO ) Lower, ( CBLAS_TRANSPOSE ) NoTrans, 
						     tile_size_mm, TILE_SIZE, 
						     -1.0,  TILE_A[mmi][ki], TILE_SIZE, 
						      1.0, TILE_A[mmi][mmi], tile_size_mm );
#else
				cblas_dsyrk( CblasRowMajor, 
					         ( CBLAS_UPLO ) Lower, ( CBLAS_TRANSPOSE ) NoTrans, 
						     tile_size_mm, TILE_SIZE, 
						     -1.0,  TILE_A[mmi][ki], TILE_SIZE, 
						      1.0, TILE_A[mmi][mmi], tile_size_mm );
#endif
			
				for ( ni = ki + 1; ni < mmi; ni++ )
				{
				
					#pragma oss task in( TILE_A[mmi][ki] ) \
                		in( TILE_A[ni][ki] ) \
                    	inout( TILE_A[mmi][ni] ) \
						shared( TILE_A ) \
 						firstprivate( ni, mmi, ki )
#ifdef LASs_WITH_CHAMELEON
						interop_dgemm(( CBLAS_TRANSPOSE ) NoTrans,
								     ( CBLAS_TRANSPOSE ) Trans,
								     TILE_SIZE, 
								     TILE_SIZE, 
								     TILE_SIZE,
							 	     -1.0, TILE_A[mmi][ki], TILE_SIZE,
							 			   TILE_A[ni][ki],  TILE_SIZE,
							 	      1.0, TILE_A[mmi][ni], TILE_SIZE );
#else
					    cblas_dgemm( CblasRowMajor, 
								     ( CBLAS_TRANSPOSE ) NoTrans,
								     ( CBLAS_TRANSPOSE ) Trans,
								     TILE_SIZE, 
								     TILE_SIZE, 
								     TILE_SIZE,
							 	     -1.0, TILE_A[mmi][ki], TILE_SIZE,
							 			   TILE_A[ni][ki],  TILE_SIZE,
							 	      1.0, TILE_A[mmi][ni], TILE_SIZE );
#endif
				}
			}
		}
	}
	else
	{
		// --UPLO = Upper--	
		for ( ki = 0; ki < nt; ki++ )
		{
			tile_size_k = ddss_tile_size( N, ki );
			
			#pragma oss task inout( TILE_A[ki][ki] ) \
				shared( TILE_A ) \
				firstprivate( ki )
#ifdef LASs_WITH_CHAMELEON
			interop_dpotrf(UPLO,
        					 tile_size_k, TILE_A[ki][ki], 
							 tile_size_k );
#else
			LAPACKE_dpotrf_work( LAPACK_ROW_MAJOR, 
							 'U',
        					 tile_size_k, TILE_A[ki][ki], 
							 tile_size_k );
#endif

			for ( mi = ki + 1; mi < nt; mi++ )
			{
				tile_size_m = ddss_tile_size( N, mi );
			
				#pragma oss task in( TILE_A[ki][ki] ) \
					inout(TILE_A[ki][mi]) \
					shared( TILE_A ) \
					firstprivate( mi, ki )
#ifdef LASs_WITH_CHAMELEON
				interop_dtrsm(( CBLAS_SIDE ) Left, ( CBLAS_UPLO ) Upper, 
					     	 ( CBLAS_TRANSPOSE ) Trans, ( CBLAS_DIAG ) NonUnit, 
						 	 TILE_SIZE, tile_size_m, 
						 	 1.0, TILE_A[ki][ki], TILE_SIZE, 
							  	  TILE_A[ki][mi], tile_size_m );
#else
				cblas_dtrsm( CblasRowMajor, 
					     	 ( CBLAS_SIDE ) Left, ( CBLAS_UPLO ) Upper, 
					     	 ( CBLAS_TRANSPOSE ) Trans, ( CBLAS_DIAG ) NonUnit, 
						 	 TILE_SIZE, tile_size_m, 
						 	 1.0, TILE_A[ki][ki], TILE_SIZE, 
							  	  TILE_A[ki][mi], tile_size_m );
#endif
			}
			for ( mmi = ki + 1; mmi < nt; mmi++ )
			{
				tile_size_mm = ddss_tile_size( N, mmi );
			
				#pragma oss task in( TILE_A[ki][mmi] ) \
					inout( TILE_A[mmi][mmi] ) \
					shared( TILE_A ) \
					firstprivate( mmi, ki )
#ifdef LASs_WITH_CHAMELEON
				interop_dsyrk(( CBLAS_UPLO ) Upper, ( CBLAS_TRANSPOSE ) Trans, 
						     tile_size_mm, TILE_SIZE, 
						     -1.0,  TILE_A[ki][mmi], tile_size_mm, 
						      1.0, TILE_A[mmi][mmi], tile_size_mm );
#else
				cblas_dsyrk( CblasRowMajor, 
					         ( CBLAS_UPLO ) Upper, ( CBLAS_TRANSPOSE ) Trans, 
						     tile_size_mm, TILE_SIZE, 
						     -1.0,  TILE_A[ki][mmi], tile_size_mm, 
						      1.0, TILE_A[mmi][mmi], tile_size_mm );
#endif
			
				for ( ni = ki + 1; ni < mmi; ni++ )
				{
					#pragma oss task in( TILE_A[ki][ni] ) \
                		in( TILE_A[ki][mmi] ) \
                    	inout( TILE_A[ni][mmi] ) \
						shared( TILE_A ) \
 						firstprivate( ni, mmi, ki )
#ifdef LASs_WITH_CHAMELEON
						interop_dgemm(( CBLAS_TRANSPOSE ) Trans,
								     ( CBLAS_TRANSPOSE ) NoTrans,
								     TILE_SIZE, 
								     TILE_SIZE, 
								     TILE_SIZE,
							 	     -1.0, TILE_A[ki][ni],  TILE_SIZE,
							 		       TILE_A[ki][mmi], TILE_SIZE,
							 	      1.0, TILE_A[ni][mmi], TILE_SIZE );
#else
					    cblas_dgemm( CblasRowMajor, 
								     ( CBLAS_TRANSPOSE ) Trans,
								     ( CBLAS_TRANSPOSE ) NoTrans,
								     TILE_SIZE, 
								     TILE_SIZE, 
								     TILE_SIZE,
							 	     -1.0, TILE_A[ki][ni],  TILE_SIZE,
							 		       TILE_A[ki][mmi], TILE_SIZE,
							 	      1.0, TILE_A[ni][mmi], TILE_SIZE );
#endif
				}
			}
		}
	}

	// --From tile data layout to flat data layout--	
	ddss_dsymtiled2flat( N, N, A, LDA, nt, nt, TILE_A, UPLO );

	// --Tile A matrix free--
	free( TILE_A );

	return Success;

}
