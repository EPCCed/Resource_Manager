#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <malloc.h>
#include <mkl.h>
#include "TreeFunctions.h"

struct Mtx_Obj 
{
	int struct_vec_index;
	int son_index;
	int i;
	int j;
	int m;
	int n;
	int lda;
	int parent;
	int isNull;
	int isLeaf;
	union
	{
		int v_values_index;
		struct Mtx_Obj * sons;
	};
};

void write_matrix_canonical( char *name, int m, int n, double *A, int lda, int v_values_index, int row, int col, FILE *fp );

void wfile_matrix_canvalues( char *name, int total_blocks, double *values, struct Mtx_Obj *Matrix_blocks, int nrows, char *file_name);

void SetCoordinates(int total_blocks, struct Mtx_Obj * Matrix_blocks, int n_levels, int *hgthtab, int *levels_sizes);

void Add_norm_to_diagonal(struct Mtx_Obj *Matrix_blocks, double *values, double norm, int total_blocks);

void Read_matrix(char *file_name, struct Mtx_Obj **Matrix_blocks, int *total_blocks, int *n_levels, int **levels_sizes, int **treetab, int **chldtab, int **rchldtab, int **nmchtab,  int **brthtab,  int **rbrthtab, int **hgthtab, int **null_bls_read, double **values, int *norm_to_add );
