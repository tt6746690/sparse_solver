#ifndef __GRAPH_H__
#define __GRAPH_H__

#include <vector>

// depth first search on sparse graph `G`
//
//  Inputs:
//  G :     A square sparse matrix
//              Note (u,v) \in G <=> G_{v,u} > 0  (! not same as adj matrix)
//  B :     A sparse column vector where nonzeros are the starting nodes
//  
//  Outputs:
//  reachset : nodes reachable in `G` from any nonzero indices in `b`
template <typename T>
void reachset(
    const CSC<T>& G,
    const CSC<T>& B,
    std::vector<int>& reachset);
template <typename T>
void reachset(
    int n,
    const int* Gp, const int* Gi, const T* Gx,
    const int* Bp, const int* Bi, const T* Bx,
    std::vector<int>& reachset);



// Implementations

#include <stack>
#include <cassert>
#include <algorithm>

#include "minitrace.h"

template <typename T>
void reachset(
    const CSC<T>& G,
    const CSC<T>& B,
    std::vector<int>& reachset)
{
    assert(B.n == 1);
    assert(G.m == G.n);
    ::reachset(G.n, G.p, G.i, G.x, B.p, B.i, B.x, reachset);
}

template <typename T>
void reachset(
    int n,
    const int* Gp, const int* Gi, const T* Gx,
    const int* Bp, const int* Bi, const T* Bx,
    std::vector<int>& reachset)
{
    MTR_SCOPE_FUNC();

    int v, p, top = 0;
    auto S = new int[n];
    auto D = std::vector<bool>(n, false);

    for (p = Bp[0]; p < Bp[1]; ++p)
        S[top++] = Bi[p];

    MTR_BEGIN("reachset", "while");
    while (top != 0) {
        v = S[--top];
        if (!D[v]) {
            D[v] = true;
            for (p = Gp[v]; p < Gp[v+1]; ++p) {
                if (!D[Gi[p]])
                    S[top++] = Gi[p];
            }
        }
    }
    MTR_END("reachset", "while");

    MTR_BEGIN("reachset", "sort");
    for (int i = 0; i < D.size(); ++i) {
        if (D[i])
            reachset.push_back(i);
    }
    MTR_END("reachset", "sort");

    delete[] S;
}


#endif 