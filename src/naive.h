#ifndef __NAIVE_H__
#define __NAIVE_H__
#include <cstddef>

// Lower triangular solver Lx=b
//
//     Inputs:
//     n: matrix dimension
//     Lp: column pointer
//     Li: row index
//     Lx: nonzeros
//     x: sparse rhs
int lsolve(int n, int* Lp, int* Li, double* Lx, double *x);

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