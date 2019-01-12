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

int repetitions = 10;
vector<int> n_threads = { 1 };
const char* trace_output = "build/trace.json";
string matr = "small";
string solver = "all";
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
            auto n_threads_str = string(argv[3]);
            auto found = n_threads_str.find('-');
            if (found != string::npos) {
                auto start = stoi(n_threads_str.substr(0, found));
                auto end   = stoi(n_threads_str.substr(found+1, n_threads_str.size()));
                n_threads.clear();
                for (int i = start; i <= end; ++i) n_threads.push_back(i);
            } else {
                n_threads[0] = stoi(n_threads_str);
            }
        }
        case 3: {
            assert(lsolve_types.find(string(argv[2])) != lsolve_types.end() || string(argv[2]) == "all");
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



    // omp_set_num_threads(n_threads);

    for (auto&& n_thread: n_threads) {
        omp_set_num_threads(n_thread);
        for (int i = 0; i < repetitions; ++i) {
            if (solver == "all") {
                for (auto&& type : lsolve_types) {
                    lsolve(type.second, L, b, x, n_thread);
                }
            } else {
                lsolve(lsolve_types.at(solver), L, b, x, n_thread);
            }
        }
    }

    // for easier plotting ...
    printf("lsolve_simple, symbolic,        1, %.16f\n", 0.);
    printf("lsolve_eigen, symbolic,         1, %.16f\n", 0.);
    printf("lsolve_eigen_par, symbolic,     1, %.16f\n", 0.);

    printf("num_threads: %d\n", omp_get_max_threads());
    mtr_shutdown();
    return 0;
}

