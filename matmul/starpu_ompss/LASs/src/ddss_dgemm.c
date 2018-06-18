#include "../include/ddss.h"

/**
 *
 * 	@file ddss_dgemm.c
 *
 * 	@brief LASs-DDSs ddss_dgemm routine.
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

/**
 *  
 *	@ingroup DDSS
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
 *			and the number of rows of the matrix C. 
 *			M must be greater than zero.
 *
 *	@param[in]
 *	N       int.
 *          N specifies the number of columns of the matrix B 
 *			and the number of columns of the matrix C.
 *		    N must be greater than zero.
 *    
 *	@param[in]
 *	K		int.
 *	       	K specifies the number of columns of the matrix A 
 *			and the number of rows of the matrix B.
 *		   	K must be greater than zero.
 *   
 *	@param[in]
 *	ALPHA   double.
 *   		ALPHA specifies the scalar alpha.
 *   
 *	@param[in]
 *	A  		double *.
 *          A is a pointer to a matrix of dimension Ma ( rows ) by Ka  
 *			( columns ), where Ma is M and Ka is K when TRANS_A = NoTrans, 
 *			and Ma is K and Ka is M otherwise.
 *
 *  @param[in]
 *	LDA     int.
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
 *	LDB     int.
 * 		    LDB specifies the number of columns of B ( row-major order ).
 *          When TRANS_B = NoTrans then LDB must be at least max( 1, N ), 
 *			otherwise LDB must be at least max( 1, K ).
 *
 *	@param[in]
 *	BETA    double.
 *   
 *	@param[in,out]
 *	C  		double *.
 *	    	C is a pointer to a matrix of dimension M by N.
 *          On exit, C is overwritten by the M by N 
 *			matrix ( ALPHA*op( A )*op( B ) + BETA*C ).
 *   
 *	@param[in]
 *	LDC     int.
 * 	    	LDC specifies the number of columns of C ( row-major order ).
 * 	    	LDC must be at least max( 1, N ).
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
 *	@sa kdgemm
 *
 **/ 

int ddss_dgemm( enum DDSS_TRANS TRANS_A, enum DDSS_TRANS TRANS_B,
	    int M, int N, int K,
        double ALPHA, double *A, int LDA,
                      double *B, int LDB,
        double BETA,  double *C, int LDC )
{
	
	// Local variables
	int An, Bn;  
	
	// Argument checking
	if ( ( TRANS_A != NoTrans ) && ( TRANS_A != Trans ) )
	{
		fprintf( stderr, "Illegal value of TRANS_A, in ddss_dgemm code\n" );
		return NoSuccess;
	}

	if ( ( TRANS_B != NoTrans ) && ( TRANS_B != Trans ) )
	{
		fprintf( stderr, "Illegal value of TRANS_B, in ddss_dgemm code\n" );
		return NoSuccess;
	}

	if ( M < 0 )
	{
		fprintf( stderr, "Illegal value of M, in ddss_dgemm code\n" );
		return NoSuccess;
	}
	
	if ( N < 0 )
	{
		fprintf( stderr, "Illegal value of N, in ddss_dgemm code\n" );
		return NoSuccess;
	}
	
	if ( K < 0 )
	{
		fprintf( stderr, "Illegal value of K, in ddss_dgemm code\n" );
		return NoSuccess;
	}
	
	if ( TRANS_A == NoTrans )
	{
		An = K;
	}	
	else
	{
		An = M;
	}

	if ( LDA < MAX( 1, An ) )
	{
		fprintf( stderr, "Illegal value of LDA, in ddss_dgemm code\n" );	
		return NoSuccess;
	}

	if ( TRANS_B == NoTrans )
	{
		Bn = N;
	}	
	else
	{
		Bn = K;
	}

	if ( LDB < MAX( 1, Bn ) )
	{
		fprintf( stderr, "Illegal value of LDB, in ddss_dgemm code\n" );	
		return NoSuccess;
	}

	if ( LDC < MAX( 1, N ) )
	{
		fprintf( stderr, "Illegal value of LDC, in ddss_dgemm code\n" );	
		return NoSuccess;
	}

	// Quick return
	if ( M == 0 || N == 0 || ( ( ALPHA == 0.0 || K == 0 ) && BETA == 1.0 ) )
	{
		return Success;
	}	

	return kdgemm( TRANS_A, TRANS_B, M, N, K, 
		   			(const double) ALPHA, A, LDA, 
                						  B, LDB, 
		   			(const double)  BETA, C, LDC );
	
}


