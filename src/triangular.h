#ifndef __TRIANGULAR_H__
#define __TRIANGULAR_H__

//  Simple implementation
//  Lower triangular solver Lx=b
//
//     Inputs:
//     n: matrix dimension
//     Lp: column pointer
//     Li: row index
//     Lx: nonzeros
//     x: sparse rhs
int lsolve_simple(int n, int* Lp, int* Li, double* Lx, double* x);


// Eigen's implementation
int lsolve_eigen(int n, int* Lp, int* Li, double* Lx, double* x);


#endif