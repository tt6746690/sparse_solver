#ifndef __TRIANGULAR_H__
#define __TRIANGULAR_H__

#include "formats.h"
 
//
//  Lower triangular solver Lx=b
//
//  Inputs:
//      t : one of lsolve_type{simple, eigen, reachset}
//      L : sparse lower triangular matrix
//      b : sparse rhs
//  Outputs:
//      x : solution to `Lx = b`
enum class lsolve_type {
    simple,
    eigen,
    reachset,
};
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

    switch (type) {
        case lsolve_type::simple:
            lsolve_simple(L.n, L.p, L.i, L.x, xvec.data());
            break;
        case lsolve_type::eigen:
            lsolve_eigen(L.n, L.p, L.i, L.x, xvec.data());
            break;
        case lsolve_type::reachset:
            lsolve_reachset(L, b, xvec.data());
            break;
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

    std::vector<int> reachset;
    ::reachset(L, B, reachset);

    int i, j, p;
    for (i = 0; i < reachset.size(); ++i) {
        j = reachset[i];
        x[j] /= Lx[Lp[j]];
        for (p = Lp[j]+1; p < Lp[j+1]; p++) {
            x[Li[p]] -= Lx[p] * x[j];
        }
    }
}




#endif