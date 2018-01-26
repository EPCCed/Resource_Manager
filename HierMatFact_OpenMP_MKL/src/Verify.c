#include "Verify.h"

void Read_cc_A_matrix_data(double *matrix, char *matrix_file_name, int m)
{
	FILE *fp;
	char buffer[200];
	int counter, i, j, k;

	// Open file where matrix is stored
	fp = fopen(matrix_file_name, "rb");
 	if (fp == NULL)
	{
	   printf ("ERROR opening given matrix data given file!");
	}

	// Fill matrix
	while(fgets(buffer, 200, fp) == NULL) {}
	i = atoi(buffer);
	for(k=2; (k<=m*m*3) && (fgets(buffer, 200, fp) != NULL); k++)
	{
		if(k%3==0) // A matrix value
		{
			matrix[i+j*m] = atof(buffer);
		}
		else
		{
			if(k%3==1) // i coord
			{
				i = atoi(buffer);
			}
			else // j coord
			{
				j = atoi(buffer);
			}
		}
	}
	if (fclose(fp)!=0) printf("Couldn't close given matrix data file.\n");
}

void Read_cc_LU_matrix_data(double *L, double *U, char *matrix_file_name, int m)
{
	FILE *fp;
	char buffer[200];
	int counter, i, j, k;

	// Open file where matrix is stored
	fp = fopen(matrix_file_name, "rb");
 	if (fp == NULL)
	{
	   printf ("ERROR opening given matrix data given file!");
	}

	// Fill matrix
	while(fgets(buffer, 200, fp) == NULL) {}
	i = atoi(buffer);
	for(k=2; (k<=m*m*3) && (fgets(buffer, 200, fp) != NULL); k++)
	{
		if(k%3==0) // A matrix value
		{
			if(j<i) L[i+j*m] = atof(buffer);
			else U[i+j*m] = atof(buffer);
		}
		else
		{
			if(k%3==1) // i coord
			{
				i = atoi(buffer);
			}
			else // j coord
			{
				j = atoi(buffer);
			}
		}
	}
	for(i=0; i<m; i++)
	{
		for(j=0; j<i; j++)
		{
			U[i+j*m] = 0.0;
		}
		L[i+i*m] = 1.0;
		for(j=i+1; j<m; j++)
		{
			L[i+j*m] = 0.0;
		}
	}
	if (fclose(fp)!=0) printf("Couldn't close given matrix data file.\n");
}

void Verify( int m, char *A_file_name, char *LU_file_name )
{
 	double *A, *LU, *L, *U, ALPHA = 1.0, BETA = 0.0, error = 0.0;
 	char TRANSA = 'N', TRANSB = 'N';
 	int i, j;

 	// Init matrices
 	A = malloc(m*m*sizeof(double));
 	L = malloc(m*m*sizeof(double));
 	U = malloc(m*m*sizeof(double));

	// Read from files and fill A, LU matrices
 	Read_cc_A_matrix_data(A, A_file_name, m);
 	Read_cc_LU_matrix_data(L, U, LU_file_name, m);

	// Perform L·U product 
	LU = malloc(m*m*sizeof(double));   
	dgemm(&TRANSA, &TRANSB, &m, &m, &m, &ALPHA, L, &m, U, &m, &BETA, LU, &m);    

	// Get error
	for (i=0; i<m*m; i++)
	{
		error += fabs(A[i] - LU[i]);
		if ( (fabs(A[i] - LU[i])) > 0.0001 )
			printf("Error at position (%d, %d)\n", i%m, i/m);
	}

	printf("Error = %22.16g\n", error);

	// Free memory
	free(A); free(L); free(U); free(LU);
}
