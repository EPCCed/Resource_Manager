#include "HLU_OmpSs.h"

// Count the amount of left brother blocks that a given Matrix_blocks[index] block has 
int CountLeftBlocks(int index, int *brthtab)
{
	int left_blocks = 0, next_left_brth = brthtab[index];
	while( next_left_brth > -1 )
	{
		left_blocks ++;
		next_left_brth = brthtab[next_left_brth];
	}
	return left_blocks;
}

// OmpSs Task to perform LU of a block (given by its index)
#pragma oss task inout( (A)[0] ) label(LU)
void LU_Task(double *A, int lu_start, int i_block, struct Mtx_Obj * Matrix_blocks, int m_size, int lu_init)
{
	int M = (int) Matrix_blocks[i_block].m, N = (int) Matrix_blocks[i_block].n, LDA = (int) Matrix_blocks[i_block].lda, INFO;
	int *IPIV = malloc(M*M*sizeof(int));
	dgetrf(&M, &N, A, &LDA, IPIV, &INFO);
	free(IPIV);
}

// OmpSs Task to perform LU TRSM to a block (given by its index) that is in the same column of given Matrix_blocks[i_block] block after performing an LU
#pragma oss task in( (u_values)[0] ) inout( (trsm_values)[0] ) label(TRSM Col)
void TRSM_Col_Task(double *u_values, double *trsm_values, int i_trsm, int i_origin, int i_u_block, int trsm_start, int u_start, struct Mtx_Obj * Matrix_blocks, int m_u, int n_u, int m_trsm, int n_trsm, int u_size, int u_lda, int trsm_size, int trsm_lda, int trsm_init, int u_init)
{
	int M = (int) Matrix_blocks[i_origin].m, N = (int) Matrix_blocks[i_origin].n, NRHS = (int) Matrix_blocks[i_origin].m, LDA = (int) Matrix_blocks[i_u_block].lda, LDB = (int) Matrix_blocks[i_trsm].lda, INFO;
	double ALPHA = 1.0;
	char SIDE = 'R', UPLO = 'U', TRANSA = 'N', DIAG = 'N';
	dtrsm(&SIDE, &UPLO, &TRANSA, &DIAG, &M, &N, &ALPHA, u_values, &LDA, trsm_values, &LDB);
}

// OmpSs Task to perform LU TRSM to a block (given by its index) that is in the same row of given Matrix_blocks[i_block] block after performing an LU
#pragma oss task in( (l_values)[0] ) inout( (trsm_values)[0] ) label(TRSM Row)
void TRSM_Row_Task(double *l_values, double *trsm_values, int i_trsm, int i_origin, int i_l_block, int trsm_start, int l_start, struct Mtx_Obj * Matrix_blocks, int m_l, int n_l, int m_trsm, int n_trsm, int l_size, int l_lda, int trsm_size, int trsm_lda, int trsm_init, int l_init)
{
	int M = (int) Matrix_blocks[i_origin].m, N = (int) Matrix_blocks[i_origin].n, NRHS = (int) Matrix_blocks[i_origin].m, LDA = (int) Matrix_blocks[i_l_block].lda, LDB = (int) Matrix_blocks[i_trsm].lda, INFO;
	char SIDE = 'L', UPLO = 'L', TRANSA = 'N', DIAG = 'U';
	double ALPHA = 1.0;
 	dtrsm(&SIDE, &UPLO, &TRANSA, &DIAG, &M, &N, &ALPHA, l_values, &LDA, trsm_values, &LDB);
}

// OmpSs Task to UPDATE values of a block (given by its index) after performing an LU
#pragma oss task in( (l_values)[0], (u_values)[0] ) inout( (upd_values)[0] ) label(UPD)
void UPD_Task(double *l_values, double *u_values, double *upd_values, int i_l_block, int i_u_block, int i_upd_block, int i_origin_block, int l_start, int u_start, int upd_start, struct Mtx_Obj * Matrix_blocks, int m_l, int n_l, int m_u, int n_u, int m_upd, int n_upd, int l_size, int l_lda, int u_size, int u_lda, int upd_size, int upd_lda, int upd_init, int l_init, int u_init)
{
	int M = (int) Matrix_blocks[i_origin_block].m, N = (int) Matrix_blocks[i_origin_block].n, K = (int) Matrix_blocks[i_origin_block].m;
	int LDA = (int) Matrix_blocks[i_l_block].lda, LDB = (int) Matrix_blocks[i_u_block].lda, LDC = (int) Matrix_blocks[i_upd_block].lda, parent_index;
	double ALPHA = -1.0, BETA = 1.0;
	char TRANSA = 'N', TRANSB = 'N';
	dgemm(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, l_values, &LDA, u_values, &LDB, &BETA, upd_values, &LDC);    
	if((int) Matrix_blocks[i_upd_block].isNull)
	{
		Matrix_blocks[i_upd_block].isNull = 0;
	}

}

// Call as much necessary TRSM_Row tasks as brother-blocks are in the same row of Matrix_blocks[i_block]
void TRSM_Row(double *values, int i_block, int i, struct Mtx_Obj Origin_block, struct Mtx_Obj * Matrix_blocks, int *rbrthtab, int matrix_m, int *sym_null_bl)
{
	int origin_bl_size, u_index, u_start, aux_u_start, row, TRSM_index, trsm_start, l_index, l_start, upd_index, upd_start, bl, aux_upd_start, aux_trsm_start, aux_l_start, m_l, n_l, m_u, n_u, m_trsm, n_trsm, m_upd, n_upd;
	int size_l, size_u, size_upd, size_trsm, lda_l, lda_u, lda_upd, lda_trsm;
	origin_bl_size = (int) Matrix_blocks[Origin_block.struct_vec_index].m * (int) Matrix_blocks[Origin_block.struct_vec_index].n;
	l_index = (int) Matrix_blocks[i_block].sons[0].struct_vec_index;
	TRSM_index = (int) Matrix_blocks[Matrix_blocks[i_block].parent].sons[i].struct_vec_index;
	u_index = (int) Matrix_blocks[Matrix_blocks[i_block].parent].sons[i].struct_vec_index;
	upd_index = (int) Matrix_blocks[Matrix_blocks[i_block].parent].sons[i].struct_vec_index;
	while ( l_index != i_block )
	{
		while ( ! (int) Matrix_blocks[l_index].isLeaf )
		{
			l_index = (int) Matrix_blocks[l_index].sons[0].struct_vec_index;
		}
		if ( ! (int) sym_null_bl[l_index] ) // L block is not null 
		{
			if ( (int) Matrix_blocks[l_index].i > (int) Matrix_blocks[l_index].j ) // Block of triang. inferior part (not diagonal)
			{
				if ( origin_bl_size != ((int) Matrix_blocks[l_index].m * (int) Matrix_blocks[l_index].n) )
				{
					l_start = 0;
					u_start = (int) Matrix_blocks[u_index].v_values_index + ((int) Matrix_blocks[l_index].j - (int) Matrix_blocks[i_block].j);
					upd_start = (int) Matrix_blocks[upd_index].v_values_index + ( (int) Matrix_blocks[l_index].i - (int) Matrix_blocks[i_block].i );
					row = 0;
					while ( l_start < ((int) Matrix_blocks[l_index].m * (int) Matrix_blocks[l_index].n) ) // Forced partitioned sub-blocks of origin block size
					{
						for ( bl = 0; bl < (int) Matrix_blocks[upd_index].m / (int) Matrix_blocks[Origin_block.struct_vec_index].m; bl++ )
						{
							aux_upd_start = upd_start + bl * (int) Matrix_blocks[upd_index].m * (int) Matrix_blocks[Origin_block.struct_vec_index].m;
							aux_u_start = u_start + bl * (int) Matrix_blocks[upd_index].m * (int) Matrix_blocks[Origin_block.struct_vec_index].m;
							m_l = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
							n_l = (int) Matrix_blocks[Origin_block.struct_vec_index].n;
							m_u = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
							n_u = (int) Matrix_blocks[Origin_block.struct_vec_index].n;
							m_upd = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
							n_upd = (int) Matrix_blocks[Origin_block.struct_vec_index].n;
							size_l = (int) Matrix_blocks[l_index].m;		lda_l = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
							size_u = (int) Matrix_blocks[u_index].m;		lda_u = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
							size_upd = (int) Matrix_blocks[upd_index].m;	lda_upd = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
							UPD_Task(&values[(int) Matrix_blocks[l_index].v_values_index + l_start], &values[aux_u_start], &values[aux_upd_start], l_index, u_index, upd_index, (int) Origin_block.struct_vec_index, (int) Matrix_blocks[l_index].v_values_index + l_start, aux_u_start, aux_upd_start, Matrix_blocks, m_l, n_l, m_u, n_u, m_upd, n_upd, size_l, lda_l, size_u, lda_u, size_upd, lda_upd, aux_upd_start, (int) Matrix_blocks[l_index].v_values_index + l_start, aux_u_start);
						}	
						row ++;
						upd_start += (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						l_start += (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						if ( row == (int) Matrix_blocks[l_index].m / (int) Matrix_blocks[Origin_block.struct_vec_index].m )
						{
							l_start += ( (int) Matrix_blocks[Origin_block.struct_vec_index].m - 1 ) * (int) Matrix_blocks[l_index].m;
							u_start += (int) Matrix_blocks[Origin_block.struct_vec_index].m;
							upd_start -= row * (int) Matrix_blocks[Origin_block.struct_vec_index].m;
							row = 0;
						}
					}
				}
				else
				{
					u_start = (int) Matrix_blocks[u_index].v_values_index + ( (int) Matrix_blocks[l_index].j - (int) Matrix_blocks[i_block].j);
					upd_start = (int) Matrix_blocks[upd_index].v_values_index + ( (int) Matrix_blocks[l_index].i - (int) Matrix_blocks[i_block].i);
					for ( bl = 0; bl < (int) Matrix_blocks[upd_index].m / (int) Matrix_blocks[Origin_block.struct_vec_index].m; bl++ )
					{
						aux_upd_start = upd_start + bl * (int) Matrix_blocks[upd_index].m * (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						aux_u_start = u_start + bl * (int) Matrix_blocks[upd_index].m * (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						m_l = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						n_l = (int) Matrix_blocks[Origin_block.struct_vec_index].n;
						m_u = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						n_u = (int) Matrix_blocks[Origin_block.struct_vec_index].n;
						m_upd = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						n_upd = (int) Matrix_blocks[Origin_block.struct_vec_index].n;
						size_l = (int) Matrix_blocks[l_index].m;		lda_l = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						size_u = (int) Matrix_blocks[u_index].m;		lda_u = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						size_upd = (int) Matrix_blocks[upd_index].m;	lda_upd = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						UPD_Task(&values[(int) Matrix_blocks[l_index].v_values_index], &values[aux_u_start], &values[aux_upd_start], l_index, u_index, upd_index, (int) Origin_block.struct_vec_index, (int) Matrix_blocks[l_index].v_values_index, aux_u_start, aux_upd_start, Matrix_blocks, m_l, n_l, m_u, n_u, m_upd, n_upd, size_l, lda_l, size_u, lda_u, size_upd, lda_upd, aux_upd_start, (int) Matrix_blocks[l_index].v_values_index, aux_u_start);
					}	
				}
			}
			else
			{
				if ( (int) Matrix_blocks[l_index].i == (int) Matrix_blocks[l_index].j ) // Diagonal block
				{
					trsm_start = (int) Matrix_blocks[TRSM_index].v_values_index + ((int) Matrix_blocks[l_index].i - (int) Matrix_blocks[i_block].i );
					for ( i = 0; i < (int) Matrix_blocks[TRSM_index].m / (int) Matrix_blocks[Origin_block.struct_vec_index].m; i++ )
					{
						aux_trsm_start = trsm_start + i * (int) Matrix_blocks[TRSM_index].m * (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						m_l = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						n_l = (int) Matrix_blocks[Origin_block.struct_vec_index].n;
						m_trsm = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						n_trsm = (int) Matrix_blocks[Origin_block.struct_vec_index].n;
						size_l = (int) Matrix_blocks[l_index].m;		lda_l = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						size_trsm = (int) Matrix_blocks[TRSM_index].m;	lda_trsm = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						TRSM_Row_Task(&values[(int) Matrix_blocks[l_index].v_values_index], &values[aux_trsm_start], TRSM_index, Origin_block.struct_vec_index, l_index, aux_trsm_start, (int) Matrix_blocks[l_index].v_values_index, Matrix_blocks, m_l, n_l, m_trsm, n_trsm, size_l, lda_l, size_trsm, lda_trsm, aux_trsm_start, (int) Matrix_blocks[l_index].v_values_index);
					}
				}
			}
		}
		if ( rbrthtab[l_index] > 0 ) // This block has brothers on the right
		{
			l_index = (int) rbrthtab[l_index]; // Next brother on the right
		}
		else
		{
			while ( (rbrthtab[l_index] < 0) && (l_index != i_block) ) // Not brothers on the right, not end 
			{
				l_index = (int) Matrix_blocks[l_index].parent;
			}
			if ( (l_index != i_block) && (rbrthtab[l_index] > 0) ) // This block has brothers on the right
			{
				l_index = (int) rbrthtab[l_index]; // Next brother on the right
			}
		}
	}
}

// Call as TRSM_Col tasks as brother-blocks are in the same column of Matrix_blocks[i_block]
void TRSM_Col(double *values, int i_block, int i, struct Mtx_Obj Origin_block, struct Mtx_Obj * Matrix_blocks, int *rbrthtab, int matrix_m, int *sym_null_bl)
{
	int origin_bl_size, u_index, u_start, row, TRSM_index, trsm_start, l_index, l_start, upd_index, upd_start, bl, aux_upd_start, aux_trsm_start, aux_l_start, m_l, n_l, m_u, n_u, m_upd, n_upd, m_trsm, n_trsm;
	int size_l, size_u, size_upd, size_trsm, lda_l, lda_u, lda_upd, lda_trsm;
	origin_bl_size = (int) Matrix_blocks[Origin_block.struct_vec_index].m * (int) Matrix_blocks[Origin_block.struct_vec_index].n;
	u_index = (int) Matrix_blocks[i_block].sons[0].struct_vec_index;
	TRSM_index = (int) Matrix_blocks[Matrix_blocks[i_block].parent].sons[i].struct_vec_index;
	l_index = (int) Matrix_blocks[Matrix_blocks[i_block].parent].sons[i].struct_vec_index;
	upd_index = (int) Matrix_blocks[Matrix_blocks[i_block].parent].sons[i].struct_vec_index;
	while ( u_index != i_block )
	{
		while ( ! (int) Matrix_blocks[u_index].isLeaf )
		{
			u_index = (int) Matrix_blocks[u_index].sons[0].struct_vec_index;
		}
		if ( ! (int) sym_null_bl[u_index] ) // U block is not null
		{
			if ( (int) Matrix_blocks[u_index].i < (int) Matrix_blocks[u_index].j ) // Block of triang. superior part (not diagonal)
			{
				if ( origin_bl_size != ((int) Matrix_blocks[u_index].m * (int) Matrix_blocks[u_index].n) )
				{
					u_start = 0;
					l_start = (int) Matrix_blocks[l_index].v_values_index + ((int) Matrix_blocks[u_index].i - (int) Matrix_blocks[i_block].i) * (int) Matrix_blocks[l_index].m;
					upd_start = (int) Matrix_blocks[upd_index].v_values_index + ( (int) Matrix_blocks[u_index].j - (int) Matrix_blocks[i_block].j ) * (int) Matrix_blocks[upd_index].m;
					row = 0;
					while ( u_start < ((int) Matrix_blocks[u_index].m * (int) Matrix_blocks[u_index].n) ) // Forced partitioned sub-blocks of origin block size
					{
						for ( bl = 0; bl < (int) Matrix_blocks[upd_index].m / (int) Matrix_blocks[Origin_block.struct_vec_index].m; bl++ )
						{
							aux_upd_start = upd_start + bl * (int) Matrix_blocks[Origin_block.struct_vec_index].m;
							aux_l_start = l_start + bl * (int) Matrix_blocks[Origin_block.struct_vec_index].m;
							m_l = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
							n_l = (int) Matrix_blocks[Origin_block.struct_vec_index].n;
							m_u = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
							n_u = (int) Matrix_blocks[Origin_block.struct_vec_index].n;
							m_upd = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
							n_upd = (int) Matrix_blocks[Origin_block.struct_vec_index].n;
							size_l = (int) Matrix_blocks[l_index].m;		lda_l = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
							size_u = (int) Matrix_blocks[u_index].m;		lda_u = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
							size_upd = (int) Matrix_blocks[upd_index].m;	lda_upd = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
							UPD_Task(&values[aux_l_start], &values[(int) Matrix_blocks[u_index].v_values_index + u_start], &values[aux_upd_start], l_index, u_index, upd_index, (int) Origin_block.struct_vec_index, aux_l_start, (int) Matrix_blocks[u_index].v_values_index + u_start, aux_upd_start, Matrix_blocks, m_l, n_l, m_u, n_u, m_upd, n_upd, size_l, lda_l, size_u, lda_u, size_upd, lda_upd, aux_upd_start, aux_l_start, (int) Matrix_blocks[u_index].v_values_index + u_start);
						}
						row ++;
						l_start += (int) Matrix_blocks[l_index].m * (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						u_start += (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						if ( row == (int) Matrix_blocks[u_index].m / (int) Matrix_blocks[Origin_block.struct_vec_index].m )
						{
							u_start += Matrix_blocks[u_index].m * (Matrix_blocks[Origin_block.struct_vec_index].m - 1);
							l_start -= (row) * (int) Matrix_blocks[l_index].m * (int) Matrix_blocks[Origin_block.struct_vec_index].m;
							upd_start += (int) Matrix_blocks[upd_index].m * (int) Matrix_blocks[Origin_block.struct_vec_index].m;
							row = 0;
						}
					}
				}
				else
				{
					l_start = (int) Matrix_blocks[l_index].v_values_index + (int) (Matrix_blocks[u_index].i - (int) Matrix_blocks[i_block].i) * (int) Matrix_blocks[l_index].m;
					upd_start = (int) Matrix_blocks[upd_index].v_values_index + ( (int) Matrix_blocks[u_index].j - (int) Matrix_blocks[i_block].j) * (int) Matrix_blocks[upd_index].m;
					for ( bl = 0; bl < (int) Matrix_blocks[upd_index].m / (int) Matrix_blocks[Origin_block.struct_vec_index].m; bl++ )
					{
						aux_upd_start = upd_start + bl * (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						aux_l_start = l_start + bl * (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						m_l = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						n_l = (int) Matrix_blocks[Origin_block.struct_vec_index].n;
						m_u = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						n_u = (int) Matrix_blocks[Origin_block.struct_vec_index].n;
						m_upd = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						n_upd = (int) Matrix_blocks[Origin_block.struct_vec_index].n;
						size_l = (int) Matrix_blocks[l_index].m;		lda_l = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						size_u = (int) Matrix_blocks[u_index].m;		lda_u = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						size_upd = (int) Matrix_blocks[upd_index].m;	lda_upd = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						UPD_Task(&values[aux_l_start], &values[(int) Matrix_blocks[u_index].v_values_index], &values[aux_upd_start], l_index, u_index, upd_index, (int) Origin_block.struct_vec_index, aux_l_start, (int) Matrix_blocks[u_index].v_values_index, aux_upd_start, Matrix_blocks, m_l, n_l, m_u, n_u, m_upd, n_upd, size_l, lda_l, size_u, lda_u, size_upd, lda_upd, aux_upd_start, aux_l_start, (int) Matrix_blocks[u_index].v_values_index);
					}	
				}
			}
			else
			{
				if ( (int) Matrix_blocks[u_index].i == (int) Matrix_blocks[u_index].j ) // Diagonal block
				{
					trsm_start = (int) Matrix_blocks[TRSM_index].v_values_index + ((int) Matrix_blocks[u_index].j - (int) Matrix_blocks[i_block].j) * (int) Matrix_blocks[TRSM_index].m;
					for ( i = 0; i < (int) Matrix_blocks[TRSM_index].m / (int) Matrix_blocks[Origin_block.struct_vec_index].m; i++ )
					{
						aux_trsm_start = trsm_start + i * (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						m_u = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						n_u = (int) Matrix_blocks[Origin_block.struct_vec_index].n;
						m_trsm = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						n_trsm = (int) Matrix_blocks[Origin_block.struct_vec_index].n;
						size_u = (int) Matrix_blocks[u_index].m;		lda_u = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						size_trsm = (int) Matrix_blocks[upd_index].m;	lda_trsm = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
						TRSM_Col_Task(&values[(int) Matrix_blocks[u_index].v_values_index], &values[aux_trsm_start], TRSM_index, (int) Origin_block.struct_vec_index, u_index, aux_trsm_start, (int) Matrix_blocks[u_index].v_values_index, Matrix_blocks, m_u, n_u, m_trsm, n_trsm, size_u, lda_u, size_trsm, lda_trsm, aux_trsm_start, (int) Matrix_blocks[u_index].v_values_index);
					}	
				}
			}
		}
		if ( rbrthtab[u_index] > 0 ) // This block has brothers on the right
		{
			u_index = (int) rbrthtab[u_index]; // Next brother on the right
		}
		else
		{
			while ( (rbrthtab[u_index] < 0) && (u_index != i_block) ) // Not brothers on the right, not end 
			{
				u_index = (int) Matrix_blocks[u_index].parent;
			}
			if ( (u_index != i_block) && (rbrthtab[u_index] > 0) ) // This block has brothers on the right
			{
				u_index = (int) rbrthtab[u_index]; // Next brother on the right
			}
		}
	}
}

// Call as much necessary UPDATE tasks as brother-blocks are under the row of Matrix_blocks[i_block] and on the right of its column
void UPD(double *values, int i_block, int i, struct Mtx_Obj Origin_block, int n_remainig_bl, struct Mtx_Obj * Matrix_blocks, int i_block_size, int *rbrthtab, int matrix_m)
{
	int origin_bl_size, UPD_index, upd_idx, upd_start, L_Block, L_index, aux_l_start, l_start, U_Block, U_index, aux_u_start, bl, aux_upd_start, k, m_l, n_l, m_u, n_u, m_upd, n_upd;
	int size_l, size_u, size_upd, lda_l, lda_u, lda_upd;
	origin_bl_size = (int) Matrix_blocks[Origin_block.struct_vec_index].m * (int) Matrix_blocks[Origin_block.struct_vec_index].n;
	UPD_index = (int) Matrix_blocks[Matrix_blocks[i_block].parent].sons[i].struct_vec_index;
	L_Block = (int) Matrix_blocks[Matrix_blocks[i_block].parent].sons[ i - ( ( (int) Matrix_blocks[UPD_index].j - (int) Matrix_blocks[i_block].j ) / i_block_size ) * ((int) Matrix_blocks[Matrix_blocks[i_block].parent].lda) ].struct_vec_index;
	U_Block = (int) Matrix_blocks[Matrix_blocks[i_block].parent].sons[ i - ( ( (int) Matrix_blocks[UPD_index].i - (int) Matrix_blocks[i_block].i ) / i_block_size ) ].struct_vec_index;
	upd_idx = UPD_index;
	while ( upd_idx <= UPD_index )
	{
		while ( ! (int) Matrix_blocks[upd_idx].isLeaf )
		{
			upd_idx = (int) Matrix_blocks[upd_idx].sons[0].struct_vec_index;
		}
		L_index = (int) Matrix_blocks[L_Block].v_values_index + (int) Matrix_blocks[upd_idx].i - (int) Matrix_blocks[UPD_index].i;
		U_index = (int) Matrix_blocks[U_Block].v_values_index + ( (int) Matrix_blocks[upd_idx].j - (int) Matrix_blocks[UPD_index].j ) * i_block_size;
		for ( bl = 0; bl < ((int) Matrix_blocks[upd_idx].m / (int) Matrix_blocks[Origin_block.struct_vec_index].m)*((int) Matrix_blocks[upd_idx].m / (int) Matrix_blocks[Origin_block.struct_vec_index].m); bl++ )
		{ // Each forced partitioned smaller size block in Matrix_blocks[upd_idx] block (if Matrix_blocks[upd_idx] block has smaller size, then will loop only 1 time)
			aux_upd_start = (int) Matrix_blocks[upd_idx].v_values_index + (int) (bl / ( (int) Matrix_blocks[upd_idx].m / (int) Matrix_blocks[Origin_block.struct_vec_index].m ) ) * Matrix_blocks[upd_idx].m * (int) Matrix_blocks[Origin_block.struct_vec_index].m +  (int) (bl % ( (int) Matrix_blocks[upd_idx].m / (int) Matrix_blocks[Origin_block.struct_vec_index].m ) ) * (int) Matrix_blocks[Origin_block.struct_vec_index].m;
			aux_u_start = U_index + (int) (bl / ((int) Matrix_blocks[upd_idx].m / (int) Matrix_blocks[Origin_block.struct_vec_index].m)) * i_block_size * (int) Matrix_blocks[Origin_block.struct_vec_index].m;
			aux_l_start = L_index + (int) (bl % ((int) Matrix_blocks[upd_idx].m / (int) Matrix_blocks[Origin_block.struct_vec_index].m)) * (int) Matrix_blocks[Origin_block.struct_vec_index].m;
			for ( k = 0; k < i_block_size / Matrix_blocks[Origin_block.struct_vec_index].m; k++ )
			{ // Each gemm that has to be performed between the L - row blocks and U - col blocks on block Matrix_blocks[upd_idx]
				m_l = (int) Matrix_blocks[Origin_block.struct_vec_index].m;	
				n_l = (int) Matrix_blocks[Origin_block.struct_vec_index].n;				
				m_u = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
				n_u = (int) Matrix_blocks[Origin_block.struct_vec_index].n;
				m_upd = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
				n_upd = (int) Matrix_blocks[Origin_block.struct_vec_index].n;
				size_l = (int) Matrix_blocks[L_Block].m;		lda_l = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
				size_u = (int) Matrix_blocks[U_Block].m;		lda_u = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
				size_upd = (int) Matrix_blocks[upd_idx].m;		lda_upd = (int) Matrix_blocks[Origin_block.struct_vec_index].m;
				UPD_Task(&values[aux_l_start + k * Matrix_blocks[Origin_block.struct_vec_index].m * i_block_size], &values[aux_u_start + k * Matrix_blocks[Origin_block.struct_vec_index].m], &values[aux_upd_start], L_Block, U_Block, upd_idx, Origin_block.struct_vec_index, aux_l_start + k * Matrix_blocks[Origin_block.struct_vec_index].m * i_block_size, aux_u_start + k * Matrix_blocks[Origin_block.struct_vec_index].m, aux_upd_start, Matrix_blocks, m_l, n_l, m_u, n_u, m_upd, n_upd, size_l, lda_l, size_u, lda_u, size_upd, lda_upd, aux_upd_start, aux_l_start + k * Matrix_blocks[Origin_block.struct_vec_index].m * i_block_size, aux_u_start + k * Matrix_blocks[Origin_block.struct_vec_index].m);
			}
		}	
		if ( rbrthtab[upd_idx] > 0 ) // This block has brothers on the right
		{
			upd_idx = (int) rbrthtab[upd_idx]; // Next brother on the right
		}
		else
		{
			while ( (rbrthtab[upd_idx] < 0) && (upd_idx <= UPD_index) )
			{  // Not brothers on the right && not end --> Try to update parent block
				upd_idx = (int) Matrix_blocks[upd_idx].parent;
			}
			if ( (upd_idx <= UPD_index) && (rbrthtab[upd_idx] > 0) )
			{ // This block has brothers on the right && not end --> Try to update next brother-block on the right
				upd_idx = (int) rbrthtab[upd_idx];
			}
		}
	} 
}

// Perform symbolic LU to guess which initial null blocks will not be null after computing the LU
void SymbolicLU(int n_blocks, int n_levels, int *nmchtab, int *brthtab, int *rbrthtab, int *hgthtab, int *levels_sizes, int **sym_null_bl, struct Mtx_Obj * Matrix_blocks, double *flops )
{
	int i_block, parent_size, n_remainig_bl, i, i_block_size, j, UPD_index, upd_index, L_Block_i, U_Block_i, son_i, total_sons, k, son_index, upd_idx, L_index, U_index, TRSM_index, bl, aux_upd_start, aux_l_start, aux_u_start, u_index, l_index;
	int size_l, size_u, size_upd, lda_l, lda_u, lda_upd, size_trsm, lda_trsm;
	int origin_bl_size, upd_start, L_Block, l_start, U_Block, m_l, n_l, m_u, n_u, m_upd, n_upd, trsm_start, u_start, row, col, m_trsm, n_trsm, aux_trsm_start;
	struct Mtx_Obj Origin_block;
	double fsize = 1.0 * (int) Matrix_blocks[0].m * (int) Matrix_blocks[0].m * (int) Matrix_blocks[0].m, fsize2;			
	for ( i_block = 0; i_block < n_blocks; i_block++ )
	{
		parent_size = (int) Matrix_blocks[Matrix_blocks[i_block].parent].m;
		
		if ( nmchtab[i_block] == 0 )
		{ // Leaf block (equivalent to use Matrix_blocks[i_block].isLeaf)
			
			if ( brthtab[i_block] == -1 || rbrthtab[i_block] == -1 || CountLeftBlocks(i_block, brthtab) % (parent_size + 1) == 0 )
			{ // 1st or last brother or diagonal block
				
				// LU
				*flops = *flops + 2.0 / 3.0 * fsize;

				// TRSM(s) and UPD(s) after LU
				if ( rbrthtab[i_block] > -1 )
			    { // Brothers on the right --> TRSM(s) + UPD(s)
					n_remainig_bl = parent_size - 1 - ((int) Matrix_blocks[i_block].son_index % parent_size);

					// TRSM(s) to brothers on the same column
			    	for( i = (int) Matrix_blocks[i_block].son_index + 1; i <= (int) Matrix_blocks[i_block].son_index + n_remainig_bl; i++ )
					{
						if ( ! (int) (*sym_null_bl)[Matrix_blocks[Matrix_blocks[i_block].parent].sons[i].struct_vec_index] ) // Block which TRSM is being applied is not null
						{
							*flops = *flops + fsize;
						}
					}

					// TRSM(s) to brothers on the same row
					for( i = (int) Matrix_blocks[i_block].son_index + ((int) Matrix_blocks[i_block].son_index % parent_size) + n_remainig_bl + 1; i <= (int) Matrix_blocks[i_block].son_index + n_remainig_bl * parent_size; i += parent_size )
			    	{
			    		//printf("TRSM Row: Matrix_blocks[i_block].son_index = %d, parent_size = %d, n_remainig_bl = %d, i = %d, \n", Matrix_blocks[i_block].son_index, parent_size, n_remainig_bl, i);
			    		if (! (int) (*sym_null_bl)[Matrix_blocks[Matrix_blocks[i_block].parent].sons[i].struct_vec_index] ) // Block which TRSM is being applied is not null	
			    		{
							*flops = *flops + fsize;
		    			}
		    		}

			    	// UPD(s) to brothers above the row and on the right of the column
					i_block_size = (int) Matrix_blocks[i_block].m;
					for( i = (int) Matrix_blocks[i_block].son_index + ((int) Matrix_blocks[i_block].son_index % parent_size) + n_remainig_bl + 2; i <= (int) Matrix_blocks[i_block].son_index + n_remainig_bl * parent_size + n_remainig_bl; i += parent_size )
			    	{
			    		for ( j = 0; j < n_remainig_bl; j++ )
			    		{
			    			UPD_index = (int) Matrix_blocks[Matrix_blocks[i_block].parent].sons[i+j].struct_vec_index;
			    			L_Block_i = (int) Matrix_blocks[Matrix_blocks[i_block].parent].sons[ i + j - ( ( (int) Matrix_blocks[UPD_index].j - (int) Matrix_blocks[i_block].j ) / i_block_size ) * ((int) Matrix_blocks[Matrix_blocks[i_block].parent].lda) ].struct_vec_index;
							U_Block_i = (int) Matrix_blocks[Matrix_blocks[i_block].parent].sons[ i + j - ( ( (int) Matrix_blocks[UPD_index].i - (int) Matrix_blocks[i_block].i ) / i_block_size ) ].struct_vec_index;
							if ( ( (!(int) (*sym_null_bl)[L_Block_i]) && (! (int) (*sym_null_bl)[U_Block_i]) ) ) // UPD block is null but (L & U not null)
							{
								(*sym_null_bl)[UPD_index] = 0;
								*flops = *flops + 2.0 * fsize;
							}
						}
			    	}
			    }
			}
		}
		else
		{ // Not a leaf block

			if ( ( brthtab[i_block] == -1 || CountLeftBlocks(i_block, brthtab) % (parent_size + 1) == 0 ) && (rbrthtab[i_block] > -1) ) // (1st or last brother or diagonal block) && (has brothers on the right)
			{ // 1st brother or diagonal block with brothers on the right --> TRSM(s) + UPD(s)
				
				fsize2 = 1.0 * levels_sizes[n_levels-hgthtab[i_block]] * levels_sizes[n_levels-hgthtab[i_block]] * levels_sizes[n_levels-hgthtab[i_block]];
				
				// For each brother on the right...
	    		n_remainig_bl = parent_size - 1 - (Matrix_blocks[i_block].son_index % parent_size);
		    	Origin_block = Matrix_blocks[i_block];
    			while ( ! Origin_block.isLeaf )
				{ // Keeps looping until Origin_block is a leaf block
					Origin_block = Origin_block.sons[0];
				}
		    	
				// TRSM(s) to brothers on the same column
		    	for( i = Matrix_blocks[i_block].son_index + 1; i <= Matrix_blocks[i_block].son_index + n_remainig_bl; i++ )
		    	{
		    		if ( ! (int) (*sym_null_bl)[Matrix_blocks[Matrix_blocks[i_block].parent].sons[i].struct_vec_index] ) // Block which TRSM is being applied is not null
					{
						*flops = *flops + fsize2;
					}
				}

				// TRSM(s) to brothers on the same row
				for( i = Matrix_blocks[i_block].son_index + (Matrix_blocks[i_block].son_index % parent_size) + n_remainig_bl + 1; i <= Matrix_blocks[i_block].son_index + (Matrix_blocks[i_block].son_index % parent_size) + n_remainig_bl + 1 + (n_remainig_bl-1) * parent_size; i+=parent_size )
		    	{
	    			if ( ! (int) (*sym_null_bl)[Matrix_blocks[Matrix_blocks[i_block].parent].sons[i].struct_vec_index] ) // Block which TRSM is being applied is not null
					{	
						*flops = *flops + fsize2;
					}
				}

		    	// UPD(s) to brothers above the row and on the right of the column
				i_block_size = (levels_sizes[n_levels-hgthtab[i_block]]);
				for( i = 0; i < n_remainig_bl; i++ )
		    	{
		    		for( j = 0; j < n_remainig_bl; j++ )
			    	{
						son_i = Matrix_blocks[i_block].son_index + (Matrix_blocks[i_block].son_index % parent_size) + n_remainig_bl + 2 + i * Matrix_blocks[Matrix_blocks[i_block].parent].lda + j;
			    		UPD_index = (int) Matrix_blocks[Matrix_blocks[i_block].parent].sons[son_i].struct_vec_index;
			    		L_Block_i = (int) Matrix_blocks[Matrix_blocks[i_block].parent].sons[ son_i - ( ( (int) Matrix_blocks[UPD_index].j - (int) Matrix_blocks[i_block].j ) / i_block_size ) * ((int) Matrix_blocks[Matrix_blocks[i_block].parent].lda) ].struct_vec_index;
						U_Block_i = (int) Matrix_blocks[Matrix_blocks[i_block].parent].sons[ son_i - ( ( (int) Matrix_blocks[UPD_index].i - (int) Matrix_blocks[i_block].i ) / i_block_size ) ].struct_vec_index;
						if ( (! (int) (*sym_null_bl)[L_Block_i]) && (! (int) (*sym_null_bl)[U_Block_i]) )
			    		{
							*flops = *flops + 2.0 * fsize2;
			    			if ( (int) (*sym_null_bl)[UPD_index] )
			    			{
			    				(*sym_null_bl)[UPD_index] = 0;
			    			}
			    			upd_idx = UPD_index;
							while ( upd_idx <= UPD_index )
							{
								while ( ! (int) Matrix_blocks[upd_idx].isLeaf )
								{
									upd_idx = (int) Matrix_blocks[upd_idx].sons[0].struct_vec_index;
								}
								L_index = (int) Matrix_blocks[L_Block_i].v_values_index + (int) Matrix_blocks[upd_idx].i - (int) Matrix_blocks[UPD_index].i;
								U_index = (int) Matrix_blocks[U_Block_i].v_values_index + ( (int) Matrix_blocks[upd_idx].j - (int) Matrix_blocks[UPD_index].j ) * i_block_size;
								for ( bl = 0; bl < ((int) Matrix_blocks[upd_idx].m / (int) Matrix_blocks[Origin_block.struct_vec_index].m)*((int) Matrix_blocks[upd_idx].m / (int) Matrix_blocks[Origin_block.struct_vec_index].m); bl++ )
								{ // Each forced partitioned smaller size block in Matrix_blocks[upd_idx] block (if Matrix_blocks[upd_idx] block has smaller size, then will loop only 1 time)
									aux_upd_start = (int) Matrix_blocks[upd_idx].v_values_index + (int) (bl / ( (int) Matrix_blocks[upd_idx].m / (int) Matrix_blocks[Origin_block.struct_vec_index].m ) ) * Matrix_blocks[upd_idx].m * (int) Matrix_blocks[Origin_block.struct_vec_index].m +  (int) (bl % ( (int) Matrix_blocks[upd_idx].m / (int) Matrix_blocks[Origin_block.struct_vec_index].m ) ) * (int) Matrix_blocks[Origin_block.struct_vec_index].m;
									aux_u_start = U_index + (int) (bl / ((int) Matrix_blocks[upd_idx].m / (int) Matrix_blocks[Origin_block.struct_vec_index].m)) * i_block_size * (int) Matrix_blocks[Origin_block.struct_vec_index].m;
									aux_l_start = L_index + (int) (bl % ((int) Matrix_blocks[upd_idx].m / (int) Matrix_blocks[Origin_block.struct_vec_index].m)) * (int) Matrix_blocks[Origin_block.struct_vec_index].m;
									for ( k = 0; k < i_block_size / Matrix_blocks[Origin_block.struct_vec_index].m; k++ )
									{ // Each gemm that has to be performed between the L - row blocks and U - col blocks on block Matrix_blocks[upd_idx]
										if ((*sym_null_bl)[upd_idx])
										{
			    							(*sym_null_bl)[upd_idx] = 0;
			    						}
									}
								}	
								if ( rbrthtab[upd_idx] > 0 ) // This block has brothers on the right
								{
									upd_idx = (int) rbrthtab[upd_idx]; // Next brother on the right
								}
								else
								{
									while ( (rbrthtab[upd_idx] < 0) && (upd_idx <= UPD_index) )
									{  // Not brothers on the right && not end --> Try to update parent block
										upd_idx = (int) Matrix_blocks[upd_idx].parent;
									}
									if ( (upd_idx <= UPD_index) && (rbrthtab[upd_idx] > 0) )
									{ // This block has brothers on the right && not end --> Try to update next brother-block on the right
										upd_idx = (int) rbrthtab[upd_idx];
									}
								}
							} 
		    			}	
		    		}
	    		}
    		}
		}
	}
}

int main (int argc, char *argv[])
{
	// ----------------------- CHECK USAGE 
	if (argc != 2)
	{
	    printf ("Usage: HLU_OmpSs <FileOfMatrix>\n");
	    return 0;
	}

	// ----------------------- VARIABLES DECLATIONS 
	int i, j, k, l, m, n_levels, n_blocks, *levels_sizes, i_block, parent_size, parent_index, i_brth, n_left_blocks, n_remainig_bl, start, UPD_index, i_first, U_Block_i, U_block_start, L_Block_i, L_block_start, N, M, K, NRHS, LDA, LDB, LDC, INFO, ret;
	int *treetab, *chldtab, *rchldtab, *nmchtab, *brthtab, *rbrthtab, *hgthtab, *null_bl_read, *sym_null_bl, size_i_block, i_values, i_read_values, i_block_size;
	int norm_to_add, brth_start, parent_lda, son_i, son_lda, n_brth, start_i, next_start, index, smaller_bl_i = 0, total_values = 0;
	char *file_name, UPLO, TRANS, TRANSA, TRANSB, DIAG, name[10], A_file_name[15], LU_file_name[15];
    struct timeval ts, te;
	struct Mtx_Obj * Matrix_blocks;
	double ALPHA, BETA, aux, *aux_values, *values, flops = 0.0;
	struct Mtx_Obj Origin_block, TRSM_Block, L_Block, U_Block, to_UPD_block, current_bl;

	// ----------------------- GET DATA FROM INPUT
	file_name = argv[1];

	// ----------------------- READ MATRIX DATA FROM FILE
	n_blocks = 0; n_levels = 0;
	Read_matrix(file_name, &Matrix_blocks, &n_blocks, &n_levels, &levels_sizes, &treetab, &chldtab, &rchldtab, &nmchtab, &brthtab, &rbrthtab, &hgthtab, &null_bl_read, &aux_values, &norm_to_add);

	// ----------------------- PERFORM SYMBOLIC LU
	sym_null_bl = malloc(n_blocks * sizeof(int));
	for(i=0; i<n_blocks; i++)
	{
		sym_null_bl[i] = null_bl_read[i];
	}
	SymbolicLU(n_blocks, n_levels, nmchtab, brthtab, rbrthtab, hgthtab, levels_sizes, &sym_null_bl, Matrix_blocks, &flops);
	
	// Generate a bigger vector with all the read elements + space for elements obtained due to UPD's
	for (i=0; i<n_blocks; i++)
	{
		if( (int) Matrix_blocks[i].isLeaf && ! sym_null_bl[i]) // Leaf block & not null
			total_values += (int) Matrix_blocks[i].m * (int) Matrix_blocks[i].n;
	}
	values = malloc(total_values * sizeof(double));
	i_values = 0; i_read_values = 0;
	for (i=0; i<n_blocks; i++)
	{
		if( (int) Matrix_blocks[i].isLeaf )
		{
			if( ! (int) Matrix_blocks[i].isNull ) // Leaf block & not null from the begining
			{
				for ( j=0; j < (int) Matrix_blocks[i].m; j++ )
				{
					for ( k=0; k < (int) Matrix_blocks[i].n; k++ )
					{
						values[i_values] = aux_values[i_read_values];
						i_values++; i_read_values++;
					}
				}
			}
			else // Leaf block & null from the begining
			{
				if (! sym_null_bl[i]) // Not null block after LU
				{
					for ( j=0; j < (int) Matrix_blocks[i].m; j++ )
					{
						for ( k=0; k < (int) Matrix_blocks[i].n; k++ )
						{
							values[i_values] = 0.0;
							i_values++;
						}
					}
				}
			}
		}
	}
	free(aux_values);

	// Set values indexes in Mtx_Obj's
	j = 0;
	for( i = 0; i < n_blocks; i++ )
	{
		if( (int) Matrix_blocks[i].isLeaf && ! (int) sym_null_bl[i]) // Leaf and not null block
		{
			Matrix_blocks[i].v_values_index = j;
			j += Matrix_blocks[i].m * Matrix_blocks[i].n;
		}
	}
	Add_norm_to_diagonal(Matrix_blocks, values, norm_to_add, n_blocks);

	sprintf(A_file_name, "A_values.txt");
	sprintf(name, "");
	wfile_matrix_canvalues(name, n_blocks, values, Matrix_blocks, levels_sizes[n_levels-1], A_file_name);

	#pragma oss taskwait
	gettimeofday(&ts, NULL);
    #pragma oss taskwait
	
	// ----------------------- PERFORM LU
	for ( i_block = 0; i_block < n_blocks; i_block++ )
	{
		parent_size = Matrix_blocks[Matrix_blocks[i_block].parent].m;
		
		if ( nmchtab[i_block] == 0 )
		{ // Leaf block (equivalent to use Matrix_blocks[i_block].isLeaf)
			
			if ( brthtab[i_block] == -1 || rbrthtab[i_block] == -1 || CountLeftBlocks(i_block, brthtab) % (parent_size + 1) == 0 )
			{ // 1st or last brother or diagonal block
				

				// LU (and PIVOTING after LU)
				if ( (int) sym_null_bl[i_block] )
				{ // LU to a null block --> ERROR!
					printf("ERROR - LU to a null block!\n");
				}
				else
				{ // LU to a not null block
					LU_Task(&values[(int) Matrix_blocks[i_block].v_values_index], (int) Matrix_blocks[i_block].v_values_index, i_block, Matrix_blocks, (int) Matrix_blocks[i_block].m, (int) Matrix_blocks[i_block].v_values_index);	
				}
				
				// TRSM(s) and UPD(s) after LU
				if ( rbrthtab[i_block] > -1 )
			    { // Brothers on the right --> TRSM(s) + UPD(s)
					n_remainig_bl = parent_size - 1 - ((int) Matrix_blocks[i_block].son_index % parent_size);

					// TRSM(s) to brothers on the same column
			    	for( i = (int) Matrix_blocks[i_block].son_index + 1; i <= (int) Matrix_blocks[i_block].son_index + n_remainig_bl; i++ )
					{
						if ( ! (int) sym_null_bl[Matrix_blocks[Matrix_blocks[i_block].parent].sons[i].struct_vec_index] ) // Block which TRSM is being applied is not null
						{
							TRSM_Col_Task(&values[(int) Matrix_blocks[i_block].v_values_index], &values[(int) Matrix_blocks[Matrix_blocks[Matrix_blocks[i_block].parent].sons[i].struct_vec_index].v_values_index], (int) Matrix_blocks[Matrix_blocks[i_block].parent].sons[i].struct_vec_index, i_block, i_block, (int) Matrix_blocks[Matrix_blocks[Matrix_blocks[i_block].parent].sons[i].struct_vec_index].v_values_index, (int) Matrix_blocks[i_block].v_values_index, Matrix_blocks, (int) Matrix_blocks[i_block].m, (int) Matrix_blocks[i_block].n, (int) Matrix_blocks[Matrix_blocks[Matrix_blocks[i_block].parent].sons[i].struct_vec_index].m, (int) Matrix_blocks[Matrix_blocks[Matrix_blocks[i_block].parent].sons[i].struct_vec_index].n, (int) Matrix_blocks[i_block].m, (int) Matrix_blocks[i_block].m, (int) Matrix_blocks[i_block].m, (int) Matrix_blocks[i_block].m, (int) Matrix_blocks[Matrix_blocks[Matrix_blocks[i_block].parent].sons[i].struct_vec_index].v_values_index, (int) Matrix_blocks[i_block].v_values_index);
						}
					}

					// TRSM(s) to brothers on the same row
					for( i = (int) Matrix_blocks[i_block].son_index + ((int) Matrix_blocks[i_block].son_index % parent_size) + n_remainig_bl + 1; i <= (int) Matrix_blocks[i_block].son_index + n_remainig_bl * parent_size; i += parent_size )
			    	{
			    		if (! (int) sym_null_bl[Matrix_blocks[Matrix_blocks[i_block].parent].sons[i].struct_vec_index] ) // Block which TRSM is being applied is not null	
			    		{
							TRSM_Row_Task(&values[(int) Matrix_blocks[i_block].v_values_index], &values[(int) Matrix_blocks[Matrix_blocks[Matrix_blocks[i_block].parent].sons[i].struct_vec_index].v_values_index], (int) Matrix_blocks[Matrix_blocks[i_block].parent].sons[i].struct_vec_index, i_block, i_block, (int) Matrix_blocks[Matrix_blocks[Matrix_blocks[i_block].parent].sons[i].struct_vec_index].v_values_index, (int) Matrix_blocks[i_block].v_values_index, Matrix_blocks, (int) Matrix_blocks[i_block].m, (int) Matrix_blocks[i_block].n, (int) Matrix_blocks[Matrix_blocks[Matrix_blocks[i_block].parent].sons[i].struct_vec_index].m, (int) Matrix_blocks[Matrix_blocks[Matrix_blocks[i_block].parent].sons[i].struct_vec_index].n, (int) Matrix_blocks[i_block].m, (int) Matrix_blocks[i_block].m, (int) Matrix_blocks[i_block].m, (int) Matrix_blocks[i_block].m, (int) Matrix_blocks[Matrix_blocks[Matrix_blocks[i_block].parent].sons[i].struct_vec_index].v_values_index, (int) Matrix_blocks[i_block].v_values_index);
		    			}
		    		}

			    	// UPD(s) to brothers above the row and on the right of the column
					i_block_size = (int) Matrix_blocks[i_block].m;
					for( i = (int) Matrix_blocks[i_block].son_index + ((int) Matrix_blocks[i_block].son_index % parent_size) + n_remainig_bl + 2; i <= (int) Matrix_blocks[i_block].son_index + n_remainig_bl * parent_size + n_remainig_bl; i += parent_size )
			    	{
			    		for ( j = 0; j < n_remainig_bl; j++ )
			    		{
			    			UPD_index = (int) Matrix_blocks[Matrix_blocks[i_block].parent].sons[i+j].struct_vec_index;
			    			L_Block_i = (int) Matrix_blocks[Matrix_blocks[i_block].parent].sons[ i + j - ( ( (int) Matrix_blocks[UPD_index].j - (int) Matrix_blocks[i_block].j ) / i_block_size ) * ((int) Matrix_blocks[Matrix_blocks[i_block].parent].lda) ].struct_vec_index;
							U_Block_i = (int) Matrix_blocks[Matrix_blocks[i_block].parent].sons[ i + j - ( ( (int) Matrix_blocks[UPD_index].i - (int) Matrix_blocks[i_block].i ) / i_block_size ) ].struct_vec_index;
							if ( (! (int) sym_null_bl[L_Block_i]) && (! (int) sym_null_bl[U_Block_i]) )
							{ // Both L, U blocks that the UPD is performed from are NOT null
								UPD_Task(&values[(int) Matrix_blocks[L_Block_i].v_values_index], &values[(int) Matrix_blocks[U_Block_i].v_values_index], &values[(int) Matrix_blocks[UPD_index].v_values_index], L_Block_i, U_Block_i, UPD_index, i_block, Matrix_blocks[L_Block_i].v_values_index, Matrix_blocks[U_Block_i].v_values_index, Matrix_blocks[UPD_index].v_values_index, Matrix_blocks, (int) Matrix_blocks[L_Block_i].m, (int) Matrix_blocks[L_Block_i].n, (int) Matrix_blocks[U_Block_i].m, (int) Matrix_blocks[U_Block_i].n, (int) Matrix_blocks[UPD_index].m, (int) Matrix_blocks[UPD_index].n, (int) Matrix_blocks[i_block].m, (int) Matrix_blocks[i_block].m, (int) Matrix_blocks[i_block].m, (int) Matrix_blocks[i_block].m, (int) Matrix_blocks[i_block].m, (int) Matrix_blocks[i_block].m, (int) Matrix_blocks[UPD_index].v_values_index, (int) Matrix_blocks[L_Block_i].v_values_index, (int) Matrix_blocks[U_Block_i].v_values_index);
							}
						}
			    	}
			    }
			}
		}
		else
		{ // Not a leaf block

			if ( ( brthtab[i_block] == -1 || CountLeftBlocks(i_block, brthtab) % (parent_size + 1) == 0 ) && (rbrthtab[i_block] > -1) ) // (1st or last brother or diagonal block) && (has brothers on the right)
			{ // 1st brother or diagonal block with brothers on the right --> TRSM(s) + UPD(s)
				
				// For each brother on the right...
	    		n_remainig_bl = parent_size - 1 - (Matrix_blocks[i_block].son_index % parent_size);
		    	Origin_block = Matrix_blocks[i_block];
    			while ( ! Origin_block.isLeaf )
				{ // Keeps looping until Origin_block is a leaf block
					Origin_block = Origin_block.sons[0];
				}
    				
		    	// TRSM(s) to brothers on the same column
		    	for( i = Matrix_blocks[i_block].son_index + 1; i <= Matrix_blocks[i_block].son_index + n_remainig_bl; i++ )
		    	{
		    		if ( ! (int) sym_null_bl[Matrix_blocks[Matrix_blocks[i_block].parent].sons[i].struct_vec_index] ) // Block which TRSM is being applied is not null
					{
						TRSM_Col(values, i_block, i, Origin_block, Matrix_blocks, rbrthtab, levels_sizes[n_levels-1], sym_null_bl);
					}
				}

				// TRSM(s) to brothers on the same row
				for( i = Matrix_blocks[i_block].son_index + (Matrix_blocks[i_block].son_index % parent_size) + n_remainig_bl + 1; i <= Matrix_blocks[i_block].son_index + (Matrix_blocks[i_block].son_index % parent_size) + n_remainig_bl + 1 + (n_remainig_bl-1) * parent_size; i+=parent_size )
		    	{
	    			if ( ! (int) sym_null_bl[Matrix_blocks[Matrix_blocks[i_block].parent].sons[i].struct_vec_index] ) // Block which TRSM is being applied is not null
					{	
						TRSM_Row(values, i_block, i, Origin_block, Matrix_blocks, rbrthtab, levels_sizes[n_levels-1], sym_null_bl);
					}
				}
		    	
		    	// UPD(s) to brothers above the row and on the right of the column
				for( i = 0; i < n_remainig_bl; i++ )
		    	{
		    		for( j = 0; j < n_remainig_bl; j++ )
			    	{
			    		son_i = Matrix_blocks[i_block].son_index + (Matrix_blocks[i_block].son_index % parent_size) + n_remainig_bl + 2 + i * Matrix_blocks[Matrix_blocks[i_block].parent].lda + j;
			    		UPD_index = (int) Matrix_blocks[Matrix_blocks[i_block].parent].sons[son_i].struct_vec_index;
			    		L_Block_i = (int) Matrix_blocks[Matrix_blocks[i_block].parent].sons[ son_i - ( ( (int) Matrix_blocks[UPD_index].j - (int) Matrix_blocks[i_block].j ) / (levels_sizes[n_levels-hgthtab[i_block]]) ) * ((int) Matrix_blocks[Matrix_blocks[i_block].parent].lda) ].struct_vec_index;
						U_Block_i = (int) Matrix_blocks[Matrix_blocks[i_block].parent].sons[ son_i - ( ( (int) Matrix_blocks[UPD_index].i - (int) Matrix_blocks[i_block].i ) / (levels_sizes[n_levels-hgthtab[i_block]]) ) ].struct_vec_index;
			    		if ( (! (int) sym_null_bl[L_Block_i]) && (! (int) sym_null_bl[U_Block_i]) )
			    		{
			    			UPD(values, i_block, son_i, Origin_block, n_remainig_bl, Matrix_blocks, levels_sizes[n_levels-hgthtab[i_block]], rbrthtab, levels_sizes[n_levels-1]);
			    		}	
	    			}
	    		}
    		}
		}
	}
	#pragma oss taskwait
	gettimeofday(&te, NULL);

	// ----------------------- SHOW RESULT (print LU [= B] matrix and execution time)
	printf("Execution time:\n");
	printf("%g\n", (te.tv_sec - ts.tv_sec)*1000 + (te.tv_usec - ts.tv_usec)/1000.0);
	printf("GFlops:\n%g\n", (flops / ((te.tv_sec - ts.tv_sec)*1000.0 + (te.tv_usec - ts.tv_usec)/1000.0) / 1000.0 )/1000.0);

	// ----------------------- VERIFY RESULT
	sprintf(LU_file_name, "LU_values.txt");
	sprintf(name, "");
	wfile_matrix_canvalues(name, n_blocks, values, Matrix_blocks, levels_sizes[n_levels-1], LU_file_name);
	Verify(levels_sizes[n_levels-1], A_file_name, LU_file_name);

	// ----------------------- FREE MEMORY
	free(treetab); free(chldtab); free(rchldtab); free(nmchtab); free(brthtab); free(rbrthtab); free(hgthtab);
	free(sym_null_bl); free(null_bl_read);
	free(values); free(Matrix_blocks);

	// ----------------------- END OF EXECUTION
	printf("\nEnd of execution!\n");
	return 0;
}
