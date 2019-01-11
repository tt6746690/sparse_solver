#include "minitrace.h"
#include "triangular.h"

void lsolve_simple(int n, int* Lp, int* Li, double* Lx, double* x) {
    MTR_SCOPE_FUNC();
    int p = 0, j = 0;

    std::cout<<"n="<<n<<'\n';
    for (int i = 0; i < n; ++i) {
        if (x[i] != 0) {
            std::cout<<i<<": "<<x[i]<<'\n';
        }
    }
    int epsilon = 1;
    n = 21000;
    for (j = 0; j < n; j++) {


        if (x[j]/Lx[Lp[j]] > epsilon || x[j]/Lx[Lp[j]] < -epsilon) {
            std::cout<<j<<": "<<x[j]<<" / "<<Lx[Lp[j]]<<" = "<<x[j]/Lx[Lp[j]]<<"\n";
        }

        x[j] /= Lx[Lp[j]];
        
        for (p = Lp[j]+1; p < Lp[j+1]; p++) {

            if (Lx[p] * x[j] > epsilon) {
                std::cout<<"\t\tx["<<Li[p]<<"] -= L["<<Li[p]<<", "<<j<<"] * x["<<j<<"] ( = "<<x[Li[p]]<< " - "<<Lx[p]<<" * "<<x[j]<<") = "<<(x[Li[p]]-Lx[p]*x[j])<<'\n';
            }
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