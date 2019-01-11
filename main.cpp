#include <omp.h>

#include "src/minitrace.h"
#include "src/multiply.h"
#include "src/triangular.h"
#include "src/formats.h"

#include <unordered_map>
#include <cassert>
#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <string>

using namespace std;


int repetitions = 1;
unordered_map<string, lsolve_type> solvers = {
    {"simple", lsolve_type::simple},
    {"eigen", lsolve_type::eigen},
    {"reachset", lsolve_type::reachset},
    {"all", lsolve_type::simple}
};
const char* trace_output = "build/trace.json";
string matr = "small";
string solver = "simple";
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
        case 4: {
            repetitions = stoi(argv[3]);
        }
        case 3: {
            assert(solvers.find(argv[2]) != solvers.end());
            solver = argv[2];
        }
        case 2: {
            matr = argv[1];
            assert(std::find(matrices.begin(), matrices.end(), matr) != matrices.end());
        }
        case 1: break;
        default: break;
    }

    Lp = "./data/"+matr+"L.mtx";
    bp = "./data/"+matr+"b.mtx";

    L = CSC<double>(Lp.c_str(), true);
    b = CSC<double>(bp.c_str(), true);
    assert(L.m == L.n);

    for (int i = 0; i < repetitions; ++i) {
        if (solver == "all") {
            lsolve(lsolve_type::simple, L, b, x);
            lsolve(lsolve_type::eigen, L, b, x);
            lsolve(lsolve_type::reachset, L, b, x);
            lsolve(lsolve_type::eigen_par, L, b, x);
        } else {
            lsolve(solvers.at(solver), L, b, x);
        }
    }


    printf("num_threads: %d\n", omp_get_num_threads());
    mtr_shutdown();
    return 0;
}

