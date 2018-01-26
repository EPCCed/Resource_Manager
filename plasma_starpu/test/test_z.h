/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @precisions normal z -> s d c
 *
 **/
#ifndef TEST_Z_H
#define TEST_Z_H

#include "test.h"

//==============================================================================
// test routines
//==============================================================================
void test_zgemm(param_value_t param[], bool run);
void test_zgeqrf(param_value_t param[], bool run);
void test_zgetrf(param_value_t param[], bool run);
void test_zpotrf(param_value_t param[], bool run);
void test_ztrsm(param_value_t param[], bool run);

#endif // TEST_Z_H
