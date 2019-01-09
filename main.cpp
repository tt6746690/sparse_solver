#include "src/minitrace.h"
#include "src/multiply.h"
#include "src/triangular.h"
#include "src/formats.h"

#include <cassert>
#include <vector>
#include <iostream>
#include <string>
using namespace std;

const char* trace_output = "build/trace.json";
string matr;
const vector<string> matrices = {"small", "torso", "tsopf"};

vector<double> x;
string Lp, bp;

void shrink_mtxs() {
    for (const auto& matr : matrices) {
        for (const auto& type : {"L", "b"}) {
            auto from = "./data/"+matr+type+".mtx";
            auto to   = "./data/s_"+matr+type+".mtx";
            auto M = COO<double>(from.c_str(), true);
            saveMatrixMarket(to.c_str(), M);
        }
    }
}

int main(int argc, char const *argv[])
{
    mtr_init(trace_output);
    if (argc == 2) {
        matr = argv[1];
        assert(std::find(matrices.begin(), matrices.end(), matr) != matrices.end());
    } else {
        matr = "small";
    }
    
    Lp = "./data/"+matr+"L.mtx";
    bp = "./data/"+matr+"b.mtx";

    auto L = CSC<double>(Lp.c_str(), true);
    auto b = CSC<double>(bp.c_str(), true);
    assert(L.m == L.n);

    csc_to_vec(b, x);
    lsolve_simple(L.n, L.p, L.i, L.x, x.data());

    printf("solution: \n");
    for (int i = 0; i < L.n; ++i) {
        if (x[i] != 0) printf("x[%d]=%f\n", i, x[i]);
    }

    csc_to_vec(b, x);
    lsolve_eigen(L.n, L.p, L.i, L.x, x.data());

    printf("solution: \n");
    for (int i = 0; i < L.n; ++i) {
        if (x[i] != 0) printf("x[%d]=%f\n", i, x[i]);
    }

    mtr_shutdown();
    return 0;
}

