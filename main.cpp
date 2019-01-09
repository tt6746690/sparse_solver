#include "src/naive.h"
#include "src/formats.h"

#include <iostream>
#include <string>
using namespace std;

int main(int argc, char const *argv[])
{

    string torsoL = "./data/torso1/torso1.mtx";
    string torsob = "./data/torso1/b_for_torso1.mtx";

    string Lp = "./data/small/L.mtx";
    string bp = "./data/small/b.mtx";

    auto L = CSC<double>(Lp.c_str(), true);
    auto b = CSC<double>(bp.c_str(), true);
    show_csc(L);
    show_csc(b);
    assert(L.m == L.n);

    vector<double> x;
    csc_to_vec(b, x);
    lsolve(L.n, L.p, L.i, L.x, x.data());

    printf("solution: \n");
    for (int i = 0; i < d; ++i) {
        if (x[i] != 0) printf("x[%d]=%f\n", i, x[i]);
        else           printf("x[%d]=0\n", i);
    }


    return 0;
}

