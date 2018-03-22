#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <float.h>
#include <mkl.h>

#include "starpu_dgemm.h"

// Number of multiplications in GEMM
#define FMULS_GEMM(m_, n_, k_) ((m_) * (n_) * (k_))
// Number of additions in GEMM
#define FADDS_GEMM(m_, n_, k_) ((m_) * (n_) * (k_))
// Flops in DGEMM 
#define FLOPS_DGEMM(m_, n_, k_) ( FMULS_GEMM((double)(m_), (double)(n_), \
	(double)(k_)) + FADDS_GEMM((double)(m_), (double)(n_), (double)(k_)) )

//#define CHECK_RESULT

int main( int argc, char const *argv[])
{

	if(argc < 4 || argc > 4)
	{ 		
		fprintf(stderr, "Usage: ./test_dgemm M N K\n" );
		return -1;
	}

	// Local variables
	int i, j;
	int m, n, k;
	double alpha; 
	double beta;
	double error;
	double max_error;
	double count_error;	
	double *A;
	double *B;
	double *C;
	double *C_test;
	struct timeval start, end;
	double flops;
	double time_ddss;
	double flops_ddss;
	double flops_ref;
	int ret;

	m = atoi( argv[1] );
	n = atoi( argv[2] );
	k = atoi( argv[3] );

	// Set seed 
	srand(time(NULL));

	max_error = 1.0;
	count_error = 0.0;

	// Checking inputs
	if ( m < 0 )
	{
		fprintf(stderr, "Illegal value of M, M must be >= 0\n");
		return -1;
	}
	if ( n < 0 )
	{
		fprintf(stderr, "Illegal value of N, N must be >= 0\n");
		return -1;
	}
	if ( k < 0 )
	{
		fprintf(stderr, "Illegal value of K, K must be >= 0\n");
		return -1;
	}

	// Matrices allocation
	A = ( double * ) malloc( sizeof( double ) * m * k );
	B = ( double * ) malloc( sizeof( double ) * k * n );
	C = ( double * ) malloc( sizeof( double ) * m * n );
#ifdef CHECK_RESULT
	C_test = ( double * ) malloc( sizeof( double ) * m * n );
#endif

	// Alpha and beta initialization
	alpha = ( double ) rand() / (double) rand() + DBL_MIN;
	beta  = ( double ) rand() / (double) rand() + DBL_MIN;
 
	// Matrix A, B, C and C_test initialization
	for ( i = 0; i < m; i++ )
	{
		for ( j = 0; j < k; j++ )
		{
			A[ i * k + j ] = ( double ) rand() / (double) rand() 
							  + DBL_MIN;
		}
	}
	
	for ( i = 0; i < k; i++ )
	{
		for ( j = 0; j < n; j++ )
		{
			B[ i * n + j ] = ( double ) rand() / (double) rand() 
							  + DBL_MIN;
		}
	}
	
	for ( i = 0; i < m; i++ )
	{
		for ( j = 0; j < n; j++ )
		{
			C[ i * n + j ] = 0.0;
#ifdef CHECK_RESULT
			C_test[ i * n + j ] = 0.0;
#endif
		}
	}
	
	interop_starpu_init();
	
	// DGEMM FLOPS
	flops = FLOPS_DGEMM( m, n, k ); 
	
	gettimeofday( &start, NULL );
	
	ret = interop_starpu_dgemm(CblasNoTrans, CblasNoTrans,
					m, n, k,
		            alpha, A, k,
		                   B, n,
		            beta, C, n );
	
	gettimeofday( &end, NULL );
	if (ret != 0) {
		return ret;
	}
		
	interop_starpu_finalize();
	
	// FLOPS achieved by the ddss_dgemm routine
	time_ddss = (double) ((end.tv_sec * 1e6 + end.tv_usec)
                                        - (start.tv_sec * 1e6 + start.tv_usec));
	flops_ddss = flops / (double) (time_ddss / 1e6);

#ifdef CHECK_RESULT
	gettimeofday( &start, NULL );

	cblas_dgemm( CblasColMajor,
				CblasNoTrans, CblasNoTrans,
				m, n, k,
				alpha, A, m,
				B, k,
				beta, C_test, m );
	
	gettimeofday( &end, NULL );
	
	// FLOPS achieved by the (reference) cblas_dgemm routine
	flops_ref = (double) flops / (((end.tv_sec * 1e6 + end.tv_usec)
				 - (start.tv_sec * 1e6 + start.tv_usec)) / 1e6);
	
	// Error computation
	for ( i = 0; i < m; i++ )
	{
		for ( j = 0; j < n; j++ )
		{
			error = abs( C[ i * n + j ] - C_test[ i * n + j ] );
			if ( max_error > error )
				max_error = error;
			count_error += error;
		}
	}

	fprintf(stdout, "Max. error = %1.2f\n", max_error );
	fprintf(stdout, "Av. error = %1.2f\n", count_error / ( m * n ) );
#endif
	
	fprintf(stdout, "Interop time = %lf us\n", time_ddss);
	fprintf(stdout, "Interop DGEMM GFLOPS = %1.2f\n", flops_ddss / 1e9 );
#ifdef CHECK_RESULT
	fprintf(stdout, "REF. DGEMM GFLOPS = %1.2f\n", flops_ref / 1e9 );
#endif
	
	return 0;

}
