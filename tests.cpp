// This tells Catch to provide a main()
//      (only do this in one cpp file)
#define CATCH_CONFIG_MAIN  
#include "catch2/catch.hpp"

#define CATCH_MARGIN 0.001
#define VECTOR_APPROX_EQUAL(A, B) \
    REQUIRE(A.size() == B.size()); \
    for (int i__ = 0; i__ < A.size(); ++i__) {    \
        if (A[i__] == -0. && B[i__] == 0.) {continue;} \
        CHECK(A[i__] == Approx(B[i__]).margin(CATCH_MARGIN));  \
    }
#define VECTOR_APPROX_EQUAL_WARN(A, B) \
    REQUIRE(A.size() == B.size()); \
    for (int i__ = 0; i__ < A.size(); ++i__) {    \
        if (A[i__] - B[i__] >= CATCH_MARGIN && A[i__] - B[i__] <= CATCH_MARGIN) { \
            WARN("lhs["<<i__<<"]="<<A[i__]<<"\trhs["<<i__<<"]="<<B[i__]<<'\n');  \
        }  \
    }


#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <functional>

#include "src/triangular.h"
#include "src/multiply.h"
#include "src/formats.h"

using namespace std;

string Lp = "./data/smallL.mtx";
string bp = "./data/smallb.mtx";

TEST_CASE("formats") {

    const auto test_coo = [](
        const string& path_to_mtx,
        int m, int n, int nnz, 
        const vector<int>& r, const vector<int>&c, const vector<double>& x)
    {
        auto M = COO<double>(path_to_mtx.c_str());
        auto expectedM = COO<double>();
        expectedM.m = m;
        expectedM.n = n;
        expectedM.nnz = nnz;
        expectedM.rowidx = (int*) malloc(nnz * sizeof(int));
        expectedM.colidx = (int*) malloc(nnz * sizeof(int));
        expectedM.values = (double*) malloc(nnz * sizeof(double));
        std::copy(r.begin(), r.end(), expectedM.rowidx);
        std::copy(c.begin(), c.end(), expectedM.colidx);
        std::copy(x.begin(), x.end(), expectedM.values);
        REQUIRE(expectedM == M);
    };

    const auto test_csc = [](
        const string& path_to_mtx,
        int m, int n, int nnz,
        const vector<int>& p, const vector<int>&i, const vector<double>& x)
    {
        auto M = CSC<double>(path_to_mtx.c_str());
        auto expectedM = CSC<double>();
        expectedM.m = m;
        expectedM.n = n;
        expectedM.nnz = nnz;
        expectedM.p = (int*) malloc((n+1) * sizeof(int));
        expectedM.i = (int*) malloc(nnz * sizeof(int));
        expectedM.x = (double*) malloc(nnz * sizeof(double));
        std::copy(p.begin(), p.end(), expectedM.p);
        std::copy(i.begin(), i.end(), expectedM.i);
        std::copy(x.begin(), x.end(), expectedM.x);
        REQUIRE(expectedM == M);
    };


    SECTION("L") {
        vector<int> p = { 0, 2, 4, 6, 9 };
        vector<int> r = { 0, 1, 1, 3, 2, 3, 1, 2, 3 };
        vector<int> c = { 0, 0, 1, 1, 2, 2, 3, 3, 3 };
        vector<double> x = { 1., 7., 2., 9., 1., 8., 4., 1., 1. };
        test_coo(Lp,4,4,9,r,c,x);
        test_csc(Lp,4,4,9,p,r,x);
    }

    SECTION("b") {
        vector<int> p = { 0, 3 };
        vector<int> r = { 0, 1, 3 };
        vector<int> c = { 0, 0, 0 };
        vector<double> x = { 0.1428571429, 3., 10. };
        test_coo(bp,4,1,3,r,c,x);
        test_csc(bp,4,1,3,p,r,x);
    }
}


const auto triangularL = [](const std::string& matr) { return "./data/"+matr+"L.mtx"; };
const auto triangularb = [](const std::string& matr) { return "./data/"+matr+"b.mtx"; };
const auto triangularsol = [](const std::string& matr) { return "./data/"+matr+"x.mtx"; };

using lsolveT = function<void (int, int*, int*, double*, double*)>;

enum class verification_method {
    backslash, mult,
};

const auto test_lsolve = [](lsolve_type type, const string& matr, verification_method meth)
{
    auto Lp = triangularL(matr);
    auto bp = triangularb(matr);
    auto L = CSC<double>(Lp.c_str(), true);
    auto b = CSC<double>(bp.c_str(), true);
    CSC<double> x;
    vector<double> xvec, yvec, bvec, solxvec;

    auto solx = CSC<double>(triangularsol(matr).c_str(), true);
    csc_to_vec(solx, solxvec);

    lsolve(type, L, b, x);

    csc_to_vec(x, xvec);
    csc_to_vec(b, bvec);

    switch (meth) {
        // 1. by verifiying with output of matlab's equivalent `x = L \ b`
        case verification_method::backslash: {
            VECTOR_APPROX_EQUAL(solxvec, xvec);
            break;
        }
        case verification_method::mult: {
            // 2. by verifying with matrix-vector multiply `y = L \cdot x` 
            //      does not work on ill-conditioned `L`
            yvec.resize(L.m, 0);
            spmv_csc(L.n, L.p, L.i, L.x, xvec.data(), yvec.data());
            VECTOR_APPROX_EQUAL(yvec, bvec);
            break;
        }
        default: break;
    }

};


TEST_CASE("triangular_small") {
    for (const auto& matr : {"s_small", "s_medium"}) {
        SECTION(matr) {
            for (auto&& type : lsolve_types) {
                test_lsolve(type.second, matr, verification_method::backslash);
                test_lsolve(type.second, matr, verification_method::mult);
            }
        }
    }
}

TEST_CASE("triangular_large", "[.][long]") {
    for (const auto& matr : {"s_torso", "s_tsopf"}) {
        SECTION("matr") {
            for (auto&& type : lsolve_types) {
                test_lsolve(type.second, matr, verification_method::backslash);
            }
        }
    }
}