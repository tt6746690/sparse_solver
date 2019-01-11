#include "multiply.h"

#include <cstdio>
#include <cstdlib>

int spmv_csc(int n, int* Ap, int* Ai, double* Ax, double* x, double* y)
{
    int p = 0, j = 0;
    if (!Ap || !x || !y) return (0) ; /* check inputs */ 

    for (j = 0; j < n; j++) {
        for (p = Ap[j]; p < Ap[j+1]; p++) {
            y[Ai[p]] += Ax[p] * x[j];
        }
    }
    return (1);
}



