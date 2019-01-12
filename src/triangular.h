#ifndef __TRIANGULAR_H__
#define __TRIANGULAR_H__

#include <unordered_map>
#include <string>
#include "utils.h"
#include "formats.h"

enum class lsolve_type {
    simple, eigen, reachset,
    eigen_par,
};
const std::unordered_map<std::string, lsolve_type> lsolve_types = {
    {"simple",      lsolve_type::simple},
    {"eigen",       lsolve_type::eigen},
    {"reachset",    lsolve_type::reachset},
    {"eigen_par",   lsolve_type::eigen_par},
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
    CSC<T>& x,
    const std::string& matr,
    int n_thread = 1);



// impl without an optimization
void lsolve_simple(int n, int* Lp, int* Li, double* Lx, double* x);

// Eigen's implementation
void lsolve_eigen(int n,int* Lp, int* Li, double* Lx, double* x);

// Solve using reachset
void lsolve_reachset_default(
    int n,  int* Lp, int* Li, double* Lx, double* x, 
    std::vector<int> reachset);
void lsolve_reachset_small(
    int n,  int* Lp, int* Li, double* Lx, double* x, 
    std::vector<int> reachset);
void lsolve_reachset_medium(
    int n,  int* Lp, int* Li, double* Lx, double* x, 
    std::vector<int> reachset);
    
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
    CSC<T>& x,
    const std::string& matr,
    int n_thread)
{
    assert(b.n == 1);
    assert(L.m == L.n);
    assert(L.n == b.m);

    std::vector<T> xvec;
    csc_to_vec(b, xvec);

    time_point_t t1, t2, t3;

    switch (type) {
        case lsolve_type::simple: {
            t1 = now();
            lsolve_simple(L.n, L.p, L.i, L.x, xvec.data());
            t2 = now();
            printf("lsolve_simple, numeric, 1, %.16f\n", diff(t1, t2));
            break;
        }
        case lsolve_type::eigen: {
            t1 = now();
            lsolve_eigen(L.n, L.p, L.i, L.x, xvec.data());
            t2 = now();
            printf("lsolve_eigen, numeric, 1, %.16f\n", diff(t1, t2));
            break;
        }
        case lsolve_type::eigen_par: {
            t1 = now();
            lsolve_eigen_par(L.n, L.p, L.i, L.x, xvec.data());
            t2 = now();
            printf("lsolve_eigen_par, numeric, %d, %.16f\n", n_thread, diff(t1, t2));
            break;
        }
        case lsolve_type::reachset: {

            t1 = now();
            std::vector<int> reachset;
            ::reachset(L, b, reachset);
            t2 = now();
            
            if (matr == "small" || matr == "s_small") {
                lsolve_reachset_small(L.n, L.p, L.i, L.x, xvec.data(), reachset);
            } else if (matr == "medium" || matr == "s_medium") {
                lsolve_reachset_medium(L.n, L.p, L.i, L.x, xvec.data(), reachset);
            } else {
                lsolve_reachset_default(L.n, L.p, L.i, L.x, xvec.data(), reachset);
            }

            t3 = now();

            printf("lsolve_reachset_%s, symbolic, 1, %.16f\n", matr.c_str(), diff(t1, t2));
            printf("lsolve_reachset_%s, numeric,  1, %.16f\n", matr.c_str(), diff(t2, t3));
            
            break;
        }
        default:
            break;
    }

    vec_to_csc(xvec, x);
}


#endif