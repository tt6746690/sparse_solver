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
    int v, p;
    auto S = std::stack<int>();
    auto D = std::vector<bool>(n, false);

    for (p = Bp[0]; p < Bp[1]; ++p)
        S.push(Bi[p]);

    while (!S.empty()) {
        v = S.top(); S.pop();
        if (!D[v]) {
            D[v] = true;
            for (p = Gp[v]; p < Gp[v+1]; ++p)
                S.push(Gi[p]);
        }
    }

    for (int i = 0; i < D.size(); ++i) {
        if (D[i])
            reachset.push_back(i);
    }
}



    




#endif 