#ifndef __TRIANGULAR_H__
#define __TRIANGULAR_H__

#include <string>
#include "utils.h"
#include "formats.h"

enum class lsolve_type {
    simple,
    eigen,
    reachset,
    eigen_par,
};
std::string lsolve_str(lsolve_type type);

//
//  Lower triangular solver Lx=b
//
//  Inputs:
//      t : one of lsolve_type{simple, eigen, reachset}
//      L : sparse lower triangular matrix
//      b : sparse rhs
//  Outputs:
//      x : solution to `Lx = b`, in dense representation
template <typename T>
void lsolve(
    lsolve_type type,
    const CSC<T>& L,
    const CSC<T>& b,
    CSC<T>& x);


// impl without an optimization
void lsolve_simple(int n, int* Lp, int* Li, double* Lx, double* x);

// Eigen's implementation
void lsolve_eigen(int n,int* Lp, int* Li, double* Lx, double* x);

// Solve using reachset
template <typename T>
void lsolve_reachset(
    const CSC<T>& L,
    const CSC<T>& b,
    T* x);
    
// simple implementation, inner loop parallelized
void lsolve_eigen_par(int n, int* Lp, int* Li, double* Lx, double* x);

// implementations 


#include <cassert>
#include <vector>
#include <iostream>

#include "minitrace.h"
#include "graph.h"


template <typename T>
void lsolve(
    lsolve_type type,
    const CSC<T>& L,
    const CSC<T>& b,
    CSC<T>& x)
{
    assert(b.n == 1);
    assert(L.m == L.n);
    assert(L.n == b.m);

    std::vector<T> xvec;
    csc_to_vec(b, xvec);

    time_point_t t1, t2;

    switch (type) {
        case lsolve_type::simple: {
            t1 = now();
            lsolve_simple(L.n, L.p, L.i, L.x, xvec.data());
            t2 = now();
            printf("lsolve_simple, %f\n", diff(t1, t2));
            break;
        }
        case lsolve_type::eigen: {
            t1 = now();
            lsolve_eigen(L.n, L.p, L.i, L.x, xvec.data());
            t2 = now();
            printf("lsolve_eigen, %f\n", diff(t1, t2));
            break;
        }
        case lsolve_type::reachset: {
            lsolve_reachset(L, b, xvec.data());
            break;
        }
        case lsolve_type::eigen_par: {
            t1 = now();
            lsolve_eigen_par(L.n, L.p, L.i, L.x, xvec.data());
            t2 = now();
            printf("lsolve_eigen_par, %f\n", diff(t1, t2));
            break;
        }
        default:
            break;
    }

    vec_to_csc(xvec, x);
}


template <typename T>
void lsolve_reachset(
    const CSC<T>& L,
    const CSC<T>& B,
    T* x)
{
    MTR_SCOPE_FUNC();

    const int* Lp;
    const int* Li;
    const int* Bp;
    const int* Bi;
    const T* Lx;
    const T* Bx;

    Lp = L.p; Li = L.i; Bp = B.p; Bi = B.i;
    Lx = L.x; Bx = B.x;

    time_point_t t1, t2, t3;

    t1 = now();
    std::vector<int> reachset;
    ::reachset(L, B, reachset);
    t2 = now();

    int i, j, p;
    for (i = 0; i < reachset.size(); ++i) {
        j = reachset[i];
        x[j] /= Lx[Lp[j]];
        for (p = Lp[j]+1; p < Lp[j+1]; p++) {
            x[Li[p]] -= Lx[p] * x[j];
        }
    }
    t3 = now();

    printf("lsolve_reachset_symbolic, %f\n", diff(t1, t2));
    printf("lsolve_reachset_numeric, %f\n", diff(t2, t3));
}


#endif