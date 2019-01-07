#ifndef __FORMAT_H__
#define __FORMAT_H__

#include "mmio.h"

// Coordinate format (COO)
//      0-based coordinate 
//      reference: https://github.com/IntelLabs/SpMP/blob/master/COO.hpp

template <typename T>
class COO {
public:
    COO();
    // `file` contains matrix market format
    COO(const char* file);
    ~COO();
public:
    // number of rows
    int m;
    // number of columns
    int n;
    // number of nonzero entries
    int nnz;
    // row indices                  size: nnz
    int* rowidx;
    // column indices               size: nnz
    int* colidx;
    // value of nonezro entries     size: nnz
    T* values;
};


// load `.mtx` -> coordinate format 
template <typename T>
void loadMatrixMarket(
    const char* file,
    COO<T>& coo);


// Compressed sparse column format (CSC)

template <typename T>
class CSC {
public:
    CSC();
    // `file` contains matrix market format
    CSC(const char* file);
    ~CSC();
public:
    // number of rows
    int m;
    // number of columns
    int n;
    // number of nonzero entries
    int nnz;
    // column pointers              size: n+1
    int* p;
    // row indices                  size: nnz
    int* i;
    // value of nonezro entries     size: nnz
    T*   x;
};

// Converts COO -> CSC
template <typename T>
void coo_to_csc(
    const COO<T>& coo,
    CSC<T>& csc);


template <typename T> 
void show_csc(const CSC<T>& csc);


// Implementation

#include <cstddef>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>

template <typename T>
CSC<T>::CSC()
    : m(0), n(0), nnz(0), p(nullptr), i(nullptr), x(nullptr)
    {};


template <typename T>
CSC<T>::CSC(const char* file) {
    auto coo = COO<T>(file);
    coo_to_csc(coo, *this);
}

template <typename T>
CSC<T>::~CSC()
{
    if (p) free(p);
    if (i) free(i);
    if (x) free(x);
}

template <typename T>
void coo_to_csc(
    const COO<T>& coo,
    CSC<T>& csc)
{
    int m = coo.m;
    int n = coo.n;
    int nnz = coo.nnz;

    csc.~CSC();
    csc.p = (int*) malloc((n+1) * sizeof(int));
    csc.i = (int*) malloc(nnz * sizeof(int));
    csc.x = (T*)   malloc(nnz * sizeof(T));

    for (int i = 0; i <= n; ++i)  csc.p[i] = 0;
    for (int i = 0; i < nnz; ++i) csc.p[coo.colidx[i]+1] += 1;
    for (int i = 0; i < n; ++i)   csc.p[i+1] += csc.p[i];

    struct ix_t {
        int i;
        T   x;
    };
    std::vector<ix_t>  ixs(nnz);

    // fill in `x` and `i`
    for (int i = 0; i < nnz; ++i) {
        int c = csc.p[coo.colidx[i]];
        ixs[c].i = coo.rowidx[i];
        ixs[c].x = coo.values[i];
        csc.p[coo.colidx[i]] += 1;
    }

    // shift back column pointers
    for (int i = n; i > 0; --i) csc.p[i] = csc.p[i-1];
    csc.p[0] = 0;

    // sort ixs
    for (int i = 0; i < n; ++i) {
        std::sort(ixs.begin()+csc.p[i], ixs.begin()+csc.p[i+1],
            [](const ix_t& a, const ix_t& b){ return a.i < b.i; });
    }

    for (int i = 0; i < nnz; ++i) {
        csc.i[i] = ixs[i].i;
        csc.x[i] = ixs[i].x;
    }

    csc.m = m;
    csc.n = n;
    csc.nnz = nnz;
}

template <typename T> 
void show_csc(const CSC<T>& csc)
{
    int n = std::min(csc.n, 10)+1;
    int nnz = std::min(csc.nnz, 10);
    printf("CSC: (%d, %d, %d)\n", csc.m, csc.n, csc.nnz);
    printf("p[:%d]: ", n);
    for (int i = 0; i < n; ++i) printf("%d ", csc.p[i]);
    printf("\n");
    printf("i[:%d]: ", nnz);
    for (int i = 0; i < nnz; ++i) printf("%d ", csc.i[i]);
    printf("\n");
    printf("x[:%d]: ", nnz);
    for (int i = 0; i < nnz; ++i) printf("%f ", (float)csc.x[i]);
    printf("\n");
}


template <typename T>
COO<T>::COO() 
    : m(0), n(0), nnz(0), rowidx(nullptr), colidx(nullptr), values(nullptr) 
    {};

template <typename T>
COO<T>::COO(const char* file) {
    loadMatrixMarket(file, *this);
}

template <typename T>
COO<T>::~COO()
{
    if (rowidx) free(rowidx);
    if (colidx) free(colidx);
    if (values) free(values);
}

template <typename T>
void loadMatrixMarket(
    const char* file,
    COO<T>& coo)
{
    FILE* fp = fopen(file, "r");
    if (fp == NULL) {
        fprintf(stderr, "Failed to open file %s\n", file);
        exit(-1);
    }

    MM_typecode matcode;
    if (mm_read_banner(fp, &matcode) != 0) {
        fprintf(stderr, "Error: could not process Matrix Market banner\n");
        exit(-1);
    }

    if (!mm_is_valid(matcode) || mm_is_dense(matcode) || !mm_is_real(matcode)) {
        fprintf(stderr, "Error: only supports real sparse matrix\n");
        exit(-1);
    }

    int m, n, nnz;
    if (mm_read_mtx_crd_size(fp, &m, &n, &nnz) != 0) {
        fprintf(stderr, "Error: could not read matrix size\n");
        exit(-1);
    }

    int* rowidx = (int*) malloc(nnz * sizeof(int));
    int* colidx = (int*) malloc(nnz * sizeof(int));
    T*   values = (T*)   malloc(nnz * sizeof(T));

    if (!values || !rowidx || !colidx) {
        fprintf(stderr, "Error: failed to allocate memory\n");
        exit(-1);
    }


    int base = 1;
    int count = 0;
    int x, y;
    double real, imag;

    while (mm_read_mtx_crd_entry(fp, &x, &y, &real, &imag, matcode) == 0) {
        // 1-based to 0-based
        x -= base;
        y -= base;

        assert(x >= 0 && x < m);
        assert(y >= 0 && y < n);

        // only care about the lower triangular 
        if (y > x) continue;

        rowidx[count] = x;
        colidx[count] = y;
        values[count] = static_cast<T>(real);

        count += 1;
    }

    assert(count <= nnz);
    fclose(fp);

    coo.m = m;
    coo.n = n;
    coo.nnz = count;
    coo.values = values;
    coo.colidx = colidx;
    coo.rowidx = rowidx;
}


#endif