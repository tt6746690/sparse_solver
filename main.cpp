#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "src/minitrace.h"
#include "src/multiply.h"
#include "src/triangular.h"
#include "src/formats.h"

#include <cassert>
#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <string>

using namespace std;
using namespace Eigen;


const char* trace_output = "build/trace.json";
string matr = "small";
const vector<string> matrices = {
    "small", "torso", "tsopf",
    "s_small", "s_medium", "s_torso", "s_tsopf"
};

vector<double> xvec, yvec, bvec;
string Lp, bp;
CSC<double> L, b, x, y;

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

    switch (argc) {
        case 2: matr = argv[1];
        case 1: break;
        default: break;
    }
    
    assert(std::find(matrices.begin(), matrices.end(), matr) != matrices.end());
    
    Lp = "./data/"+matr+"L.mtx";
    bp = "./data/"+matr+"b.mtx";

    L = CSC<double>(Lp.c_str(), true);
    b = CSC<double>(bp.c_str(), true);
    assert(L.m == L.n);

    // int rep = 5;
    // for (int i = 0; i < rep; ++i) {
    //     lsolve(lsolve_type::simple, L, b, x);
    //     lsolve(lsolve_type::eigen, L, b, x);
    //     lsolve(lsolve_type::reachset, L, b, x);
    // }


    lsolve(lsolve_type::simple, L, b, x);

    csc_to_vec(b, bvec);
    csc_to_vec(x, xvec);

    yvec.resize(L.m, 0);    
    spmv_csc(L.n, L.p, L.i, L.x, xvec.data(), yvec.data());
    vec_to_csc(yvec, y);

    for (int i = 0; i < yvec.size(); ++i) {
        if (yvec[i] != bvec[i]) {
            cout<<setprecision(numeric_limits<double>::digits10+1)
                <<i<<": (Lx=y, b) = ("<<yvec[i]<<", "<<bvec[i]<<")\n";
        }
    }

    show_csc(y);
    show_csc(b);


    mtr_shutdown();
    return 0;
}

