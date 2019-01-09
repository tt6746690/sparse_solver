// This tells Catch to provide a main()
//      (only do this in one cpp file)
#define CATCH_CONFIG_MAIN  
#include "catch2/catch.hpp"

#include <string>
#include <vector>
#include <algorithm>
#include <iostream>

#include "src/naive.h"
#include "src/formats.h"

using namespace std;

// Formats

TEST_CASE("formats") {

    const auto test_coo = [](
        const string& path_to_mtx,
        int m, int n, int nnz, 
        const vector<int>& r, const vector<int>&c, const vector<double>& x)
    {
        COO<double> M = COO<double>(path_to_mtx.c_str());
        COO<double> expectedM = COO<double>();
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
        CSC<double> M = CSC<double>(path_to_mtx.c_str());
        CSC<double> expectedM = CSC<double>();
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

    string Lp = "./data/small/L.mtx";
    string bp = "./data/small/b.mtx";

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



TEST_CASE("triangular solve") {

    string Lp = "./data/small/L.mtx";
    string bp = "./data/small/b.mtx";

    auto L = CSC<double>(Lp.c_str());
    auto b = CSC<double>(bp.c_str());

    const auto test_lsolve = [](
        const string& Lp, const string& bp,
        const vector<double> expected_x)
    {
        auto L = CSC<double>(Lp.c_str(), true);
        auto b = CSC<double>(bp.c_str(), true);
        assert(L.m == L.n);

        vector<double> x;
        csc_to_vec(b, x);
        lsolve(L.n, L.p, L.i, L.x, x.data());
        for (int i = 0; i < x.size(); ++i) {
            REQUIRE(expected_x[i] == Approx(x[i]));
        }
    };


    SECTION("small") {
        test_lsolve(Lp, bp, { 1./7., 1., 0., 1. });
    }
}