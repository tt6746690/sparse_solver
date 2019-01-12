
#include "triangular.h"
#include "minitrace.h"

void lsolve_reachset_small(
    int n,  int* Lp, int* Li, double* Lx, double* x, 
    std::vector<int> reachset)
{
MTR_SCOPE_FUNC();
        int p;
    int px;
    int j;
    for (px = 0; px < 3; ++px) {
    j = reachset[px];
    x[j] /= Lx[Lp[j]];
    for (p = Lp[j] + 1; p < Lp[j + 1]; ++p) {
    x[Li[p]] -= Lx[p] * x[j];
};
};

}