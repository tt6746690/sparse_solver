
#include "triangular.h"
#include "minitrace.h"

void lsolve_reachset_medium(
    int n,  int* Lp, int* Li, double* Lx, double* x, 
    std::vector<int> reachset)
{
MTR_SCOPE_FUNC();
        int p;
    int px;
    int j;
    x[0] /= Lx[0];
    for (p = 1; p < 3; ++p) {
    x[Li[p]] -= Lx[p] * x[0];
};
    for (px = 1; px < 3; ++px) {
    j = reachset[px];
    x[j] /= Lx[Lp[j]];
    for (p = Lp[j] + 1; p < Lp[j + 1]; ++p) {
    x[Li[p]] -= Lx[p] * x[j];
};
};
    x[7] /= Lx[20];
    for (p = 21; p < 23; ++p) {
    x[Li[p]] -= Lx[p] * x[7];
};
    for (px = 4; px < 6; ++px) {
    j = reachset[px];
    x[j] /= Lx[Lp[j]];
    for (p = Lp[j] + 1; p < Lp[j + 1]; ++p) {
    x[Li[p]] -= Lx[p] * x[j];
};
};

}