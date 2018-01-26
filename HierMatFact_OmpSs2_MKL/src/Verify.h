#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <malloc.h>
#include <mkl.h>

void Read_cc_A_matrix_data(double *matrix, char *matrix_file_name, int m);

void Read_cc_LU_matrix_data(double *L, double *U, char *matrix_file_name, int m);

void Verify( int m, char *A_file_name, char *LU_file_name );
