#ifndef __MULTIPLY_H__
#define __MULTIPLY_H__
#include <cstddef>

// Sparse matrix-vector multiply: y = A*x
//     
//     Inputs:
//     n: matrix dimension
//     Ap: column pointer
//     Ai: row index
//     Ax: nonzeros
//     x: a dense vector
//
//     Outputs:
//     y: a dense vector
int spmv_csc(int n, size_t *Ap, int *Ai, double *Ax, double *x, double *y);


#endif