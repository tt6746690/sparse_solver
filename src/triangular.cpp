#include "minitrace.h"
#include "triangular.h"

int lsolve_simple(int n, int* Lp, int* Li, double* Lx, double* x) {
    MTR_SCOPE_FUNC();
    int p = 0, j = 0;
    if (!Lp || !Li || !x) return (0) ; /* check inputs */

    for (j = 0; j < n; j++) {
        x[j] /= Lx[Lp[j]];
        for (p = Lp[j]+1; p < Lp[j+1]; p++) {
            x[Li[p]] -= Lx[p] * x[j];
        }
    }
    return (1);
}


int lsolve_eigen(int n, int* Lp, int* Li, double* Lx, double* x) {
    MTR_SCOPE_FUNC();
    int p = 0, j = 0;
    if (!Lp || !Li || !x) return (0) ; /* check inputs */

    for (j = 0; j < n; j++) {
        if (x[j] != 0) {
            x[j] /= Lx[Lp[j]];
            for (p = Lp[j]+1; p < Lp[j+1]; p++) {
                x[Li[p]] -= Lx[p] * x[j];
            }
        }
    }
    return (1);
}
