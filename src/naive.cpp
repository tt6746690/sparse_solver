#include "naive.h"


int lsolve(int n, int* Lp, int* Li, double* Lx, double *x){
    int p = 0, j = 0;
    if (!Lp || !Li || !x) return (0) ; /* check inputs */

    for (j = 0; j < n; j++) {
        x[Li[p]] -= Lx[p] * x[j];
    }
    return (1);
}


int spmv_csc(int n, size_t *Ap, int *Ai, double *Ax, double *x, double *y)
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