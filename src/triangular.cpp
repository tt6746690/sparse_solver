#include <omp.h>

#include "minitrace.h"
#include "triangular.h"

std::string lsolve_str(lsolve_type type) {
    switch (type) {
        case lsolve_type::simple:       return "simple";
        case lsolve_type::eigen:        return "eigen";
        case lsolve_type::reachset:     return "reachset";
        case lsolve_type::eigen_par:    return "eigen_par";
        default: return "";
    }
}

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
                std::cout<<"\t\tx["<<Li[p]<<"] -= L["<<Li[p]<<", "<<j<<"] * x["<<j
                <<"] ( = "<<x[Li[p]]<< " - "<<Lx[p]<<" * "<<x[j]<<") = "<<(x[Li[p]]-Lx[p]*x[j])<<'\n';
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


void lsolve_eigen_par(int n, int* Lp, int* Li, double* Lx, double* x) {
    MTR_SCOPE_FUNC();
    int p = 0, j = 0;

    for (j = 0; j < n; j++) {
        if (x[j] != 0) {
            x[j] /= Lx[Lp[j]];

#pragma omp parallel for 
            for (p = Lp[j]+1; p < Lp[j+1]; p++) {
                x[Li[p]] -= Lx[p] * x[j];
            }
        }
    }
}

void lsolve_reachset_default(
    int n,  int* Lp, int* Li, double* Lx, double* x, 
    std::vector<int> reachset)
{
    MTR_SCOPE_FUNC();
    int i, j, p;
    for (i = 0; i < reachset.size(); ++i) {
        j = reachset[i];
        x[j] /= Lx[Lp[j]];
        for (p = Lp[j]+1; p < Lp[j+1]; p++) {
            x[Li[p]] -= Lx[p] * x[j];
        }
    }
}

