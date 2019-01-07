#include "src/naive.h"
#include "src/formats.h"

#include <iostream>
#include <string>
using namespace std;

int main(int argc, char const *argv[])
{

    string torsoL = "./data/torso1/torso1.mtx";
    string torsob = "./data/torso1/b_for_torso1.mtx";

    auto L = CSC<double>(torsoL.c_str());
    auto b = CSC<double>(torsob.c_str());
    show_csc(L);
    show_csc(b);

    assert(L.m == L.n);
    int d = L.m;
    double* x = (double*) calloc(d, sizeof(double));
    for (int i = 0; i < b.nnz; ++i) x[b.i[i]] = b.x[i];

    lsolve(L.n, L.p, L.i, L.x, x);

    for (int i = 0; i < d; ++i) {
        if (x[i] != 0) 
            printf("x[%d]=%f\n", i, x[i]);
    }


    return 0;
}

