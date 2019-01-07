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
    int m;
    int n;
    int nnz;
    int* rowidx;
    int* colidx;
    T* values;
};


// load `.mtx` -> coordinate format 
template <typename T>
void loadMatrixMarket(
    const char* file,
    COO<T>& coo);


// Implementation

#include <cstddef>
#include <cstdlib>
#include <cstdio>
#include <cassert>

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

    // size of matrix

    int m, n, nnz;
    if (mm_read_mtx_crd_size(fp, &m, &n, &nnz) != 0) {
        fprintf(stderr, "Error: could not read matrix size\n");
        exit(-1);
    }

    // memory allocation

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