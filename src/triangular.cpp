#include "minitrace.h"
#include "triangular.h"

void lsolve_simple(int n, int* Lp, int* Li, double* Lx, double* x) {
    MTR_SCOPE_FUNC();
    int p = 0, j = 0;

#ifdef NUMERICAL_PRINT
    int epsilon = 1;
    n = 21000;
#endif

    for (j = 0; j < n; j++) {

#ifdef NUMERICAL_PRINT
        if (x[j]/Lx[Lp[j]] > epsilon || x[j]/Lx[Lp[j]] < -epsilon) {
            std::cout<<j<<": "<<x[j]<<" / "<<Lx[Lp[j]]<<" = "<<x[j]/Lx[Lp[j]]<<"\n";
        }
#endif

        x[j] /= Lx[Lp[j]];
        
        for (p = Lp[j]+1; p < Lp[j+1]; p++) {

#ifdef NUMERICAL_PRINT
            if (Lx[p] * x[j] > epsilon) {
                std::cout<<"\t\tx["<<Li[p]<<"] -= L["<<Li[p]<<", "<<j<<"] * x["<<j<<"] ( = "<<x[Li[p]]<< " - "<<Lx[p]<<" * "<<x[j]<<") = "<<(x[Li[p]]-Lx[p]*x[j])<<'\n';
            }
#endif

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