#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <errno.h>
#include <math.h>
#include <sys/time.h>
#include <sys/times.h>
#include <malloc.h>
#include <mkl.h>
#include "omp.h"
#include "Verify.h"
#include "Read_matrix_from_file.h" 

// Count the amount of left brother blocks that a given Matrix_blocks[index] block has 
int CountLeftBlocks(int index, int *brthtab);

// OmpSs Task to perform LU of a block (given by its index)
void LU_Task(double *A, int lu_start, int i_block, struct Mtx_Obj * Matrix_blocks, int m_size, int lu_init);

// OmpSs Task to perform LU TRSM to a block (given by its index) that is in the same column of given Matrix_blocks[i_block] block after performing an LU
void TRSM_Col_Task(double *u_values, double *trsm_values, int i_trsm, int i_origin, int i_u_block, int trsm_start, int u_start, struct Mtx_Obj * Matrix_blocks, int m_u, int n_u, int m_trsm, int n_trsm, int u_size, int u_lda, int trsm_size, int trsm_lda, int trsm_init, int u_init);

// OmpSs Task to perform LU TRSM to a block (given by its index) that is in the same row of given Matrix_blocks[i_block] block after performing an LU
void TRSM_Row_Task(double *l_values, double *trsm_values, int i_trsm, int i_origin, int i_l_block, int trsm_start, int l_start, struct Mtx_Obj * Matrix_blocks, int m_l, int n_l, int m_trsm, int n_trsm, int l_size, int l_lda, int trsm_size, int trsm_lda, int trsm_init, int l_init);

void UPD_Task(double *l_values, double *u_values, double *upd_values, int i_l_block, int i_u_block, int i_upd_block, int i_origin_block, int l_start, int u_start, int upd_start, struct Mtx_Obj * Matrix_blocks, int m_l, int n_l, int m_u, int n_u, int m_upd, int n_upd, int l_size, int l_lda, int u_size, int u_lda, int upd_size, int upd_lda, int upd_init, int l_init, int u_init);

// Call as much necessary TRSM_Row tasks as brother-blocks are in the same row of Matrix_blocks[i_block]
void TRSM_Row(double *values, int i_block, int i, struct Mtx_Obj Origin_block, struct Mtx_Obj * Matrix_blocks, int *rbrthtab, int matrix_m, int *sym_null_bl);

// Call as TRSM_Col tasks as brother-blocks are in the same column of Matrix_blocks[i_block]
void TRSM_Col(double *values, int i_block, int i, struct Mtx_Obj Origin_block, struct Mtx_Obj * Matrix_blocks, int *rbrthtab, int matrix_m, int *sym_null_bl);

// Call as much necessary UPDATE tasks as brother-blocks are under the row of Matrix_blocks[i_block] and on the right of its column
void UPD(double *values, int i_block, int i, struct Mtx_Obj Origin_block, int n_remainig_bl, struct Mtx_Obj * Matrix_blocks, int i_block_size, int *rbrthtab, int matrix_m);

// Perform symbolic LU to guess which initial null blocks will not be null after computing the LU
void SymbolicLU(int n_blocks, int n_levels, int *nmchtab, int *brthtab, int *rbrthtab, int *hgthtab, int *levels_sizes, int **sym_null_bl, struct Mtx_Obj * Matrix_blocks, double *flops );

int main (int argc, char *argv[]);
