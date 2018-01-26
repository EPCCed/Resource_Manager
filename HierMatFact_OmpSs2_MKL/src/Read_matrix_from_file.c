#include "Read_matrix_from_file.h"

void write_matrix_canonical( char *name, int m, int n, double *A, int lda, int v_values_index, int row, int col, FILE *fp )
{
  int i, j;
  for (i=0; i<m; i++)
    for (j=0; j<n; j++)
    {
    	fprintf(fp, "%d\n%d\n%22.16g\n", row+i, col+j, A[v_values_index + j*lda+i]);
	}
}

void wfile_matrix_canvalues( char *name, int total_blocks, double *values, struct Mtx_Obj *Matrix_blocks, int nrows, char *file_name)
{
	int i, j, k;
	char toWrite[100];
	FILE *fp;
	fp = fopen(file_name, "w");
	if (fp == NULL)
	{
		printf("Couldn't open file for writing data.\n");
	}
   
	for( k=0; k < total_blocks; k++ )
	{
		if ( (int) Matrix_blocks[k].isLeaf ) // Null block
		{
			if (! (int) Matrix_blocks[k].isNull ) // Not null block
				write_matrix_canonical(name, Matrix_blocks[k].m, Matrix_blocks[k].n, values, Matrix_blocks[k].lda, Matrix_blocks[k].v_values_index, Matrix_blocks[k].i, Matrix_blocks[k].j, fp);
			else
				for (i=0; i<Matrix_blocks[k].m; i++)
					for (j=0; j<Matrix_blocks[k].n; j++)
					{
						fprintf(fp, "%d\n%d\n0.0\n", Matrix_blocks[k].i+i, Matrix_blocks[k].j+j);
					}
		}
	}
	fclose(fp);	
}

// Assign (i,j) to each block
void SetCoordinates(int total_blocks, struct Mtx_Obj * Matrix_blocks, int n_levels, int *hgthtab, int *levels_sizes)
{
	int level, i;
	Matrix_blocks[total_blocks-1].i = 0;
	Matrix_blocks[total_blocks-1].j = 0;
	for(level = 1; level < n_levels; level++)
	{
		for(i=0; i<total_blocks; i++)
		{
			if(hgthtab[i]-1 == level)
			{
				Matrix_blocks[i].i = Matrix_blocks[Matrix_blocks[i].parent].i + levels_sizes[n_levels-level-1] * (int) (Matrix_blocks[i].son_index % Matrix_blocks[Matrix_blocks[i].parent].lda);
				Matrix_blocks[i].j = Matrix_blocks[Matrix_blocks[i].parent].j + levels_sizes[n_levels-level-1] * (int) (Matrix_blocks[i].son_index / Matrix_blocks[Matrix_blocks[i].parent].lda);
			}
		}
	}
}

void Add_norm_to_diagonal(struct Mtx_Obj *Matrix_blocks, double *values, double norm, int total_blocks)
{
	int i, m, j, v_index;
	for (i=0; i<total_blocks; i++)
	{
		if( ((int) Matrix_blocks[i].i == (int) Matrix_blocks[i].j) && (int) Matrix_blocks[i].isLeaf && (! (int) Matrix_blocks[i].isNull) ) // Diagonal leaf block
		{
			m = (int) Matrix_blocks[i].m;
			for (j=0; j<m*m; j+=m+1)
			{ // Diagonal elements
				v_index = (int) Matrix_blocks[i].v_values_index + j;
				values[v_index] = (double) values[v_index] + norm;
			}
		}
	}
}

void Read_matrix(char *file_name, struct Mtx_Obj **Matrix_blocks, int *total_blocks, int *n_levels, int **levels_sizes, int **treetab, int **chldtab, int **rchldtab, int **nmchtab,  int **brthtab,  int **rbrthtab, int **hgthtab, int **null_bls_read, double **values, int *norm_to_add )
{
	FILE *fp;
 	char buffer[100], *aux = "", name[10];
 	int i, j, k, *blocks_sizes, values_index = 0, added_children = 0, last_size = 0, i_val = 0, total_values = 0;

 	fp = fopen(file_name, "rb");
 	if (fp == NULL)
	{
	   printf ("ERROR opening given file!");
	}

 	// 1st element of the file: value of the norm (to add to diagonal values)
 	while(fgets(buffer, 100, fp) == NULL) {}
	*norm_to_add = atoi(buffer);

 	// 2nd element of the file: total of blocks
 	fgets(buffer, 100, fp);
	*total_blocks = atoi(buffer);
	(*Matrix_blocks) = malloc(*total_blocks * sizeof(struct Mtx_Obj));
    (*chldtab) = malloc(*total_blocks * sizeof(int));
    (*rchldtab) = malloc(*total_blocks * sizeof(int));
    (*nmchtab) = malloc(*total_blocks * sizeof(int));
    (*brthtab) = malloc(*total_blocks * sizeof(int));
    (*rbrthtab) = malloc(*total_blocks * sizeof(int));
    (*hgthtab) = malloc(*total_blocks * sizeof(int));
    (*treetab) = malloc(*total_blocks * sizeof(int));
	blocks_sizes = malloc(*total_blocks * sizeof(int));
	(*null_bls_read) = malloc( *total_blocks * sizeof(int));

	// Store every parents indexes in treetab = parents_indexes vector (stored in postorder in given file)
	for ( i=0; (i < *total_blocks) && (fgets(buffer, 100, fp) != NULL); i++ )
	{
		(*Matrix_blocks)[i].struct_vec_index = i;
		(*Matrix_blocks)[i].parent = atoi(buffer);
		(*treetab)[i] = atoi(buffer);
	}

	ComputeEliminationTreeVectors((*treetab), (*chldtab), (*rchldtab), (*nmchtab), (*brthtab), (*rbrthtab), (*hgthtab), *total_blocks);
	*n_levels = (*hgthtab)[0];
	*levels_sizes = (int *) malloc(*n_levels * sizeof(int));

	// Store every block size
	j = 0;
	for ( i=0; (i < *total_blocks) && (fgets(buffer, 100, fp) != NULL); i++ )
    {
		blocks_sizes[i] = atoi(buffer);
		if(last_size < atoi(buffer))
		{
			(*levels_sizes)[j] = sqrt(atoi(buffer));
			last_size = atoi(buffer);
			j ++;
		}
    }

   	// Store if blocks are null (1) or not (0)	
	for ( i=0; (i < *total_blocks) && (fgets(buffer, 100, fp) != NULL); i++ )
	{
		if (atoi(buffer))
		{ // Null block
			 (*Matrix_blocks)[i].isNull = 1;
		}
		else
		{
			(*Matrix_blocks)[i].isNull = 0;
		 	total_values += blocks_sizes[i];
		}
		(*null_bls_read)[i] = atoi(buffer);
	}

	// Store matrix values
	(*values) = malloc( total_values * sizeof(double));
	fgets(buffer, 100, fp);
	for ( i=0; (i < *total_blocks); i++ )
	{
		if ( (*chldtab)[i] == -1 ) // i-th block is leaf block
		{
			(*Matrix_blocks)[i].isLeaf = 1;
			(*Matrix_blocks)[i].m = sqrt(blocks_sizes[i]);
			(*Matrix_blocks)[i].n = (*Matrix_blocks)[i].m;
			(*Matrix_blocks)[i].lda = (*Matrix_blocks)[i].m;
			if ( ! (*Matrix_blocks)[i].isNull )
			{ // Not null block then store values
				for ( j=0; j < blocks_sizes[i]; j++ )
				{
					(*values)[i_val] = atof(buffer);
					fgets(buffer, 100, fp);
					i_val++;
				}
			}
		}
		else // i-th block is not a leaf block
		{
			added_children = 0;
			j = 0;
			(*Matrix_blocks)[i].isLeaf = 0;
			(*Matrix_blocks)[i].m = sqrt((*nmchtab)[i]);
			(*Matrix_blocks)[i].n = (*Matrix_blocks)[i].m;
			(*Matrix_blocks)[i].lda = (*Matrix_blocks)[i].m;
			(*Matrix_blocks)[i].sons = malloc((*nmchtab)[i] * sizeof(struct Mtx_Obj));
			while ( added_children < (*nmchtab)[i] )
			{
				if ( (*treetab)[j] == i )
				{
					(*Matrix_blocks)[i].sons[added_children] = (*Matrix_blocks)[j];
					(*Matrix_blocks)[j].son_index = added_children;
					added_children++;
				}
				j++;
			}
		}
	}

	SetCoordinates(*total_blocks, (*Matrix_blocks), *n_levels, *hgthtab, *levels_sizes);
}
