#include "../include/ddss.h"

/**
 *
 * 	@file kdgemm.c
 *
 * 	@brief LASs-DDSs kdgemm routine.
 *
 * 	LASs-DDSs is a software package provided by:
 * 	Barcelona Supercomputing Center - Centro Nacional de Supercomputacion
 *
 * 	@author Pedro Valero-Lara pedro.valero@bsc.es
 * 	@date 2017-01-05
 * 	@reviewer 
 * 	@modified 
 *
 **/

/**
 *  
 *	@ingroup DDSs
 *   
 *	Performs the matrix-matrix operation:
 *
 *		C = ALPHA * op( A ) * op( B ) + BETA * C
 *
 *	where op( X ) is one of:
 *
 *		op( X ) = X      or
 *		op( X ) = X**T 
 *	
 *	ALPHA and BETA are scalars, and A, B and C are matrices, 
 *	with op( A ) an M by K matrix, op( B ) a K by N matrix and C
 *	an M by N matrix.
 *
**/

/**
 *	
 * 	@param[in]
 *	TRANS_A	enum DDSS_TRANS.
 *	        TRANS_A specifies the form of op( A ) to be used in
 *			the matrix multiplication as follows:
 *      	- NoTrans:    op( A ) = A.
 *       	- Trans:      op( A ) = A**T.
 *	
 *	@param[in]
 *	TRANS_B	enum DDSS_TRANS.
 *	    	TRANS_B specifies the form of op( B ) to be used in
 *          the matrix multiplication as follows:
 *			- NoTrans:    op( B ) = B.
 *			- Trans:      op( B ) = B**T.
 *
 *	@param[in]
 *	M       int.
 *			M specifies the number of rows of the matrix A 
 *			and of the matrix C. 
 *			M must be greater than zero.
 *
 *	@param[in]
 *	N       int.
 *          N specifies the number of columns of the matrix B 
 *			and the number of columns of the matrix C.
 *		    N must be greater than zero.
 *    
 *	@param[in]
 *  K		int.
 *	       	K specifies the number of columns of the matrix A 
 *			and the number of rows of the matrix B.
 *		   	K must be greater than zero.
 *   
 *	@param[in]
 *  ALPHA   double.
 *   
 *	@param[in]
 *	A  		double *.
 *          A is a pointer to a matrix of dimension Ma ( rows ) by Ka  
 *			( columns ), where Ma is M and Ka is K when TRANS_A = NoTrans, 
 *			and Ma is K and Ka is M otherwise.
 *
 *  @param[in]
 *  LDA     int.
 *          LDA specifies the number of columns of A ( row-major order ). 
 *			When TRANS_A = NoTrans then LDA must be at least max( 1, K ), 
 *			otherwise LDA must be at least max( 1, M ).
 *   
 *	@param[in]
 *	B  		double *.
 *          B is a pointer to a matrix of dimension Kb ( rows ) by Nb  
 *			( columns ), where Kb is K and Nb is N when TRANS_B = NoTrans, 
 *			and Kb is N and Nb is K otherwise.
 *   
 *	@param[in]
 *  LDB     int.
 * 		    LDB specifies the number of columns of B ( row-major order ).
 *          When TRANS_B = NoTrans then LDB must be at least max( 1, N ), 
 *			otherwise LDB must be at least max( 1, K ).
 *
 *	@param[in]
 *  BETA    double.
 *   
 *	@param[in,out]
 *  C  		double *.
 *	    	C is a pointer to a matrix of dimension LDC by N.
 *          On exit, C is overwritten by the M by N 
 *			matrix ( ALPHA*op( A )*op( B ) + BETA*C ).
 *   
 *	@param[in]
 *  LDC     int.
 * 	    	LDC specifies the number of columns of C ( row-major order).
 * 	    	LDC must be at least max( 1, M )
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
 *	@sa ddss_dgemm
 *	@sa ddss_tile
 *	@sa ddss_flat2tiled
 *	@sa ddss_tiled2flat
 *
 **/ 

enum DDSS_RETURN kdgemm( enum DDSS_TRANS TRANS_A, enum DDSS_TRANS TRANS_B,
	    int M, int N, int K,
        const double ALPHA, double *A, int LDA,
                            double *B, int LDB,
        const double BETA,  double *C, int LDC )
{

	// Local variables
	int mt, kt, nt;
	int mi, ki, ni;
	int Am, An;
	int Bm, Bn;
	int tile_size_m;
	int tile_size_n;
	int tile_size_k;
	int k_check;
	double betat;
	enum DDSS_RETURN ret;

	// Number of tiles
	if ( M % TILE_SIZE == 0 )
	{
		mt = M / TILE_SIZE;
	}
	else
	{
		mt = ( M / TILE_SIZE ) + 1;
	}

	if ( K % TILE_SIZE == 0 )
	{
		kt = K / TILE_SIZE;
	}
	else
	{
		kt = ( K / TILE_SIZE ) + 1;
	}

	if ( N % TILE_SIZE == 0 )
	{
		nt = N / TILE_SIZE;
	}
	else
	{
		nt = ( N / TILE_SIZE ) + 1;
	}
	
    /****************************
	--Tile matrices declaration--
	****************************/

	if ( TRANS_A == NoTrans )
	{
		Am = mt;
		An = kt;
	}
	else
	{
		Am = kt;
		An = mt;
	}
		
	if ( TRANS_B == NoTrans )
	{
		Bm = kt;
		Bn = nt;
	}
	else
	{
		Bm = nt;
		Bn = kt;
	}

	/***************************
	--Tile matrices allocation--
	***************************/
	
	double (*TILE_A)[An][TILE_SIZE * TILE_SIZE] = malloc ( Am * An * 
		TILE_SIZE * TILE_SIZE * sizeof( double ) );

	if ( TILE_A == NULL)
	{
		fprintf( stderr, "Failure in kdgemm for matrix TILE_A\n" );
		return NoSuccess;
	}

	double (*TILE_B)[Bn][TILE_SIZE * TILE_SIZE] = malloc ( Bm * Bn * 
		TILE_SIZE * TILE_SIZE * sizeof( double ) );
    
	if ( TILE_B == NULL)
	{
		fprintf( stderr, "Failure in kdgemm for matrix TILE_B\n" );
		return NoSuccess;
	}

	double (*TILE_C)[nt][TILE_SIZE * TILE_SIZE] = malloc ( mt * nt * 
		TILE_SIZE * TILE_SIZE * sizeof( double ) );

	if ( TILE_C == NULL)
	{
		fprintf( stderr, "Failure in kdgemm for matrix TILE_C\n" );
		return NoSuccess;
	}

	/*********************************************
	--From flat data layout to tiled data layout--	
	*********************************************/

	// From flat matrix A to tile matrix TILE_A
	if ( TRANS_A == NoTrans )
	{
		ddss_dflat2tiled( M, K, A, LDA, mt, kt, TILE_A );
	}
	else
	{
		ddss_dflat2tiled( K, M, A, LDA, kt, mt, TILE_A );
	}

	// From flat matrix B to tile matrix TILE_B 
	if ( TRANS_B == NoTrans )
	{
		ddss_dflat2tiled( K, N, B, LDB, kt, nt, TILE_B );
	}
	else
	{
		ddss_dflat2tiled( N, K, B, LDB, nt, kt, TILE_B );
	}
	
	// From flat matrix C to tile matrix TILE_C 
	ddss_dflat2tiled( M, N, C, LDC, mt, nt, TILE_C );

	/*************
	--DGEMM tile--
	*************/	
	
	for ( mi = 0; mi < mt; mi++ )
	{
		tile_size_m = ddss_tile_size( M, mi );
		for ( ni = 0; ni < nt; ni++ )
		{
			tile_size_n = ddss_tile_size( N, ni );
			if ( TRANS_A == NoTrans )
			{
				k_check = K;
			}
			else
			{
				k_check = M;
			}
			// --Scale on C ( C = BETA * C )--	
			if ( ( ALPHA == 0.0 ) || ( k_check == 0 ) )
			{
				#pragma oss task inout( TILE_C[mi][ni] ) \
					shared( TILE_A, TILE_B, TILE_C ) \
					firstprivate( mi, ni ) \
					no_copy_deps
				cblas_dgemm( CblasRowMajor, 
							( CBLAS_TRANSPOSE ) TRANS_A,
							( CBLAS_TRANSPOSE ) TRANS_B,
							tile_size_m, 
							tile_size_n, 
									  0,
							ALPHA,	 TILE_A[0][0], 1, 
								     TILE_B[0][0], 1, 
							 BETA, TILE_C[mi][ni], tile_size_n );
			}
			else if ( TRANS_A == NoTrans )
			{
				// --TRANS_A = NoTrans & TRANS_B = NoTrans--	
				if ( TRANS_B == NoTrans )
				{
					for ( ki = 0; ki < kt; ki++ )
					{
						tile_size_k = ddss_tile_size( K, ki );
						if (ki == 0)
						{
							betat = BETA;
						}
						else 
						{
							betat = 1.0;
						}  
						
						#pragma oss task in( TILE_A[mi][ki] ) \
                     		in( TILE_B[ki][ni] ) \
                     		inout( TILE_C[mi][ni] ) \
							shared( TILE_A, TILE_B, TILE_C ) \
 							firstprivate( mi, ni, ki, betat ) \
							no_copy_deps
						cblas_dgemm( CblasRowMajor, 
									( CBLAS_TRANSPOSE ) TRANS_A,
									( CBLAS_TRANSPOSE ) TRANS_B,
									tile_size_m, 
									tile_size_n, 
									tile_size_k,
							 		ALPHA, TILE_A[mi][ki], tile_size_k,
							 			   TILE_B[ki][ni], tile_size_n,
							 		betat, TILE_C[mi][ni], tile_size_n );
					}
				}
				// --TRANS_A = NoTrans & TRANS_B = Trans--	
				else
				{
					for ( ki = 0; ki < kt; ki++ )
					{
						tile_size_k = ddss_tile_size( K, ki );
						if (ki == 0)
						{
							betat = BETA;
						}
						else 
						{
							betat = 1.0;
						}  
							
						#pragma oss task in( TILE_A[mi][ki] ) \
                     		in( TILE_B[ni][ki] ) \
                     		inout( TILE_C[mi][ni] ) \
							shared( TILE_A, TILE_B, TILE_C ) \
 							firstprivate( mi, ni, ki, betat )
						cblas_dgemm( CblasRowMajor, 
									( CBLAS_TRANSPOSE ) TRANS_A,
									( CBLAS_TRANSPOSE ) TRANS_B,
									tile_size_m, 
									tile_size_n, 
									tile_size_k,
							 		ALPHA,	TILE_A[mi][ki], tile_size_k,
							 				TILE_B[ni][ki], tile_size_k,
							 		betat,  TILE_C[mi][ni], tile_size_n );
					}
				}
			}
			else
			{
				// --TRANS_A = Trans & TRANS_B = NoTrans--	
				if ( TRANS_B == NoTrans )
				{
					for ( ki = 0; ki < kt; ki++ )
					{
						tile_size_k = ddss_tile_size( K, ki );
						if (ki == 0)
						{
							betat = BETA;
						}
						else 
						{
							betat = 1.0;
						}  

						#pragma oss task in( TILE_A[ki][mi] ) \
                     		in( TILE_B[ki][ni] ) \
                     		inout( TILE_C[mi][ni] ) \
							shared( TILE_A, TILE_B, TILE_C ) \
 							firstprivate( mi, ni, ki, betat )
						cblas_dgemm( CblasRowMajor, 
                                    ( CBLAS_TRANSPOSE ) TRANS_A,
									( CBLAS_TRANSPOSE ) TRANS_B,
									tile_size_m, 
									tile_size_n, 
									tile_size_k,
							 		ALPHA,	TILE_A[ki][mi], tile_size_m,
							 				TILE_B[ki][ni], tile_size_n,
							 		betat,  TILE_C[mi][ni], tile_size_n );
					}
				}	
				// --TRANS_A = Trans & TRANS_B = Trans--	
				else
				{
					for ( ki = 0; ki < kt; ki++ )
					{
						tile_size_k = ddss_tile_size( K, ki );
						if (ki == 0)
						{
							betat = BETA;
						}
						else 
						{
							betat = 1.0;
						}  
						
						#pragma oss task in( TILE_A[ki][mi] ) \
                     		in( TILE_B[ni][ki] ) \
                     		inout( TILE_C[mi][ni] ) \
							shared( TILE_A, TILE_B, TILE_C ) \
 							firstprivate( mi, ni, ki, betat )
						cblas_dgemm( CblasRowMajor, 
                                    ( CBLAS_TRANSPOSE ) TRANS_A,
									( CBLAS_TRANSPOSE ) TRANS_B,
									tile_size_m, 
									tile_size_n, 
									tile_size_k,
							 		ALPHA,	TILE_A[ki][mi], tile_size_m,
							 				TILE_B[ni][ki], tile_size_k,
							 		betat,  TILE_C[mi][ni], tile_size_n );
					}
				}
			}
		}
	}
	
	/************************************************
	// --From tiled data layout to flat data layout--	
	************************************************/
	
	// From tile matrix TILE_C to flat matrix C 
	ddss_dtiled2flat( M, N, C, LDC, mt, nt, TILE_C );

	// --Tile matrices free--
	free( TILE_A );
	free( TILE_B );
	free( TILE_C );

	return Success;

}
