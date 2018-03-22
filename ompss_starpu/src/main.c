#include "../include/ddss.h"
#include "../include/interop.h"

/**
 *
 * 	@file ddss_test_dpotrf.c
 *
 * 	@brief LASs-DDSs ddss_test_dgemm routine.
 *
 * 	LASs-DDSs is a software package provided by:
 * 	Barcelona Supercomputing Center - Centro Nacional de Supercomputacion
 *
 * 	@author Pedro Valero-Lara pedro.valero@bsc.es
 * 	@date 2017-08-16
 * 	@reviewer 
 * 	@modified 
 *
 **/

/**
 *  
 *	@ingroup DDSs
 *   
 *  main:
 *	Performs the testing of the ddss_dpotrf routine. 
 *
**/

/**
 *
 *	@sa ddss_dpotrf
 *	@sa kdpotrf
 *
 **/

//#define CHECK_RESULT

int main( int argc, char const *argv[])
{

	if( argc < 3 || argc > 3 )
	{ 		
		fprintf(stderr, "Usage: ./test_dpotrf N UPLO\n" );
		return NoSuccess;
	}

	
	// Local variables
	int i, j;
	int n;
	const char *uplo_input = NULL;
	enum DDSS_UPLO uplo = Lower;
	double error;
	double max_error;
	double count_error;	
	double *A;
	double *A_test;
	struct timeval start, end;
	double flops;
	double flops_ddss; 
	double flops_ref; 
	enum DDSS_RETURN ret;
	lapack_int retval;
	int seed[] = {0, 0, 0 , 1};

	n = atoi( argv[1] );
	
	if ( strlen( argv[2] ) != 1 ) 
	{
		fprintf(stderr,"Illegal value of UPLO, UPLO can be L or U\n");
		return NoSuccess;
	}
	uplo_input = argv[2];	
	
	max_error = 0.0;
	count_error = 0.0;

	// Checking inputs
	if ( n < 0 )
	{
		fprintf(stderr, "Illegal value of N, N must be >= 0\n");
		return NoSuccess;
	}

	if (uplo_input[0] == 'U' )
	{
		uplo = Upper;
	}
	else if ( uplo_input[0] == 'L' )
	{
		uplo = Lower;
	}
	else
	{
		fprintf(stderr, "Illegal value of UPLO, UPLO can be U or L\n");
		return NoSuccess;
	}
	
	// Matrices allocation
	A = ( double * ) malloc( sizeof( double ) * n * n );
	A_test = ( double * ) malloc( sizeof( double ) * n * n );

	// Matrix A initialization
    retval = LAPACKE_dlarnv( 1, seed, n * n, A );
    assert( retval == 0 );
	
	// Make A and A_test matrices symmetric positive definite
	for ( i = 0; i < n; i++ )
	{
		A[ i * n + i ] += n; 
		for ( j = 0; j < n; j++ )
		{
			A_test[ j * n + i ] = A[ i * n + j ];
		}
	}

	memcpy( A_test, A, n * n * sizeof( double ) );

#ifdef LASs_WITH_CHAMELEON
	interop_starpu_init();
#endif

	// DPOTRF FLOPS
	flops = FLOPS_DPOTRF( n ); 
	
	gettimeofday( &start, NULL );
	
	ret = ddss_dpotrf( uplo, n, A, n );
	
	gettimeofday( &end, NULL );
	if ( ret == NoSuccess )
		return ret;
	
#ifdef LASs_WITH_CHAMELEON
	interop_starpu_finalize();
#endif

	// FLOPS achieved by the ddss_dgemm routine
	flops_ddss = flops / (double) (((end.tv_sec * 1e6 + end.tv_usec) 
					- (start.tv_sec * 1e6 + start.tv_usec)) / 1e6);

#ifdef CHECK_RESULT
	gettimeofday( &start, NULL );

	LAPACKE_dpotrf_work( LAPACK_ROW_MAJOR, 
						 uplo_input[0],
        				 n, A_test, n );
	
	gettimeofday( &end, NULL );
	
	// FLOPS achieved by the (reference) cblas_dgemm routine
	flops_ref = (double) flops / (((end.tv_sec * 1e6 + end.tv_usec)
				 - (start.tv_sec * 1e6 + start.tv_usec)) / 1e6);
	
	// Error computation
	for ( i = 0; i < n; i++ )
	{
		for ( j = 0; j < n; j++ )
		{
			error = fabs( A[ i * n + j ] - A_test[ i * n + j ] );
			if ( max_error > error )
				max_error = error;
			count_error += error;
		}
	}

	fprintf(stdout, "Max. error = %1.2f\n", max_error );
	fprintf(stdout, "Av. error = %1.2f\n", count_error / ( n * n ) );
#endif

	fprintf(stdout, "Interop DPOTRF GFLOPS = %1.2f\n", flops_ddss / 1e9 );
#ifdef CHECK_RESULT
	fprintf(stdout, "REF. DPOTRF GFLOPS = %1.2f\n", flops_ref / 1e9 );
#endif
	
	return Success;

}
