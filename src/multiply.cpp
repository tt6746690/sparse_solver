#include "multiply.h"

#include <cstdio>
#include <cstdlib>

int spmv_csc(int n, int* Ap, int* Ai, double* Ax, double* x, double* y)
{
    int p = 0, j = 0;
    if (!Ap || !x || !y) return (0) ; /* check inputs */ 

    for (j = 0; j < n; j++) {
        // printf("%d-th column:\n p: ", j);
        for (p = Ap[j]; p < Ap[j+1]; p++) {
            // printf("y[%d] += %f * %f", Ai[p], Ax[p], x[j]);
            y[Ai[p]] += Ax[p] * x[j];
        }
        // printf("\n after %d-th column y: \n", j);
        // for (int i = 0; i < n; ++i) {
        //     printf("y[%d]=%f\n",i, y[i]);
        // }
    }
    return (1);
}



