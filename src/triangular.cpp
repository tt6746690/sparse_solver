#include "minitrace.h"
#include "triangular.h"

void lsolve_simple(int n, int* Lp, int* Li, double* Lx, double* x) {
    MTR_SCOPE_FUNC();
    int p = 0, j = 0;

    for (j = 0; j < n; j++) {
        x[j] /= Lx[Lp[j]];
        for (p = Lp[j]+1; p < Lp[j+1]; p++) {
            x[Li[p]] -= Lx[p] * x[j];
        }
    }
}


void lsolve_eigen(int n, int* Lp, int* Li, double* Lx, double* x) {
    MTR_SCOPE_FUNC();
    int p = 0, j = 0;

    for (j = 0; j < n; j++) {
        if (x[j] != 0) {
            x[j] /= Lx[Lp[j]];
            for (p = Lp[j]+1; p < Lp[j+1]; p++) {
                x[Li[p]] -= Lx[p] * x[j];
            }
        }
    }
}
