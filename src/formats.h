#ifndef __FORMAT_H__
#define __FORMAT_H__

#include "mmio.h"
#include <vector>

// Coordinate format (COO)
//      0-based coordinate 
//      reference: https://github.com/IntelLabs/SpMP/blob/master/COO.hpp

template <typename T>
class COO {
public:
    COO();
    // `file` contains matrix market format
    COO(const char* file, bool lower_triangular = false);
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

template <typename T>
inline bool operator==(const COO<T>& lhs, const COO<T>& rhs);


// load/save `.mtx` -> coordinate format 
//
//  lower_triangular: bool
//      keep nonzeros when  `rowidx <= colidx`
template <typename T>
void loadMatrixMarket(
    const char* file,
    COO<T>& coo,
    bool lower_triangular);

template <typename T>
void saveMatrixMarket(
    const char* file,
    const COO<T>& coo);


// Compressed sparse column format (CSC)

template <typename T>
class CSC {
public:
    CSC();
    // `file` contains matrix market format
    CSC(const char* file, bool lower_triangular = false);
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


template <typename T>
inline bool operator==(const CSC<T>& lhs, const CSC<T>& rhs);

// Convert sparse column vector to/from dense represeentation
template <typename T>
void csc_to_vec(
    const CSC<T>& csc,
    std::vector<T>& v);
template <typename T>
void vec_to_csc(
    const std::vector<T>& v,
    CSC<T>& csc);


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
#include <iostream>
#include <algorithm>

// COO

template <typename T>
COO<T>::COO() 
    : m(0), n(0), nnz(0), rowidx(nullptr), colidx(nullptr), values(nullptr) 
    {};

template <typename T>
COO<T>::COO(const char* file, bool lower_triangular) {
    loadMatrixMarket(file, *this, lower_triangular);
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
    COO<T>& coo,
    bool lower_triangular)
{
    FILE* fp = fopen(file, "r");
    if (fp  == NULL) {
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
        if (lower_triangular && y > x) continue;

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

template <typename T>
void saveMatrixMarket(
    const char* file,
    const COO<T>& coo)
{
    FILE *fp = fopen(file, "w");
    if (fp == NULL) {
        fprintf(stderr, "Fail to open file %s\n", file);
        exit(-1);
    }

    MM_typecode matcode;
    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_coordinate(&matcode);
    mm_set_real(&matcode);

    if (int err = mm_write_mtx_crd((char *)file, 
            coo.m, coo.n, coo.nnz, coo.rowidx, coo.colidx, coo.values, matcode, 0)) {
        fprintf(stderr, "Fail to write matrix to %s (error code = %d)\n", file, err);
    }
}

template <typename T>
inline bool operator==(const COO<T>& lhs, const COO<T>& rhs)
{
    if (lhs.m == rhs.m && lhs.n == rhs.n && lhs.nnz == rhs.nnz &&
        std::equal(lhs.rowidx, lhs.rowidx+lhs.nnz, rhs.rowidx) &&
        std::equal(lhs.colidx, lhs.colidx+lhs.nnz, rhs.colidx) &&
        std::equal(lhs.values, lhs.values+lhs.nnz, rhs.values))
    {
        return true;
    }
    return false;
}


// CSC

template <typename T>
CSC<T>::CSC()
    : m(0), n(0), nnz(0), p(nullptr), i(nullptr), x(nullptr)
    {};


template <typename T>
CSC<T>::CSC(const char* file, bool lower_triangular) {
    auto coo = COO<T>(file, lower_triangular);
    coo_to_csc(coo, *this);
}

template <typename T>
CSC<T>::~CSC()
{
#ifdef MEMORY_PRINT
    std::cout<<"CSC::~CSC()"<<": "<<this<<" ";
    if (!p && !i && !x) std::cout<<"ptrs are nullptr\n";
    else std::cout<<"releasing ptrs!\n";
#endif

    if (p) free(p);
    if (i) free(i);
    if (x) free(x);
}


template <typename T>
inline bool operator==(const CSC<T>& lhs, const CSC<T>& rhs)
{
    if (lhs.m == rhs.m && lhs.n == rhs.n && lhs.nnz == rhs.nnz &&
        std::equal(lhs.p, lhs.p+(lhs.n+1), rhs.p) &&
        std::equal(lhs.i, lhs.i+lhs.nnz, rhs.i) &&
        std::equal(lhs.x, lhs.x+lhs.nnz, rhs.x))
    {
        return true;
    }
    return false;
}

template <typename T>
void csc_to_vec(
    const CSC<T>& csc,
    std::vector<T>& v)
{
    assert(csc.n == 1);
    v.clear();
    v.resize(csc.m, 0);
    for (int i = 0; i < csc.nnz; ++i)
        v[csc.i[i]] = csc.x[i];
}

template <typename T>
void vec_to_csc(
    const std::vector<T>& v,
    CSC<T>& csc)
{
    int m = v.size();
    int n = 1;
    int nnz = 0;

    int* p = (int*) malloc((n+1) * sizeof(int));
    int* i = (int*) malloc(v.size() * sizeof(int));
    T*   x = (T*)   malloc(v.size() * sizeof(T));

    for (int j = 0; j < v.size(); ++j) {
        if (v[j] != 0) {
            i[nnz] = j;
            x[nnz] = v[j];
            nnz += 1;
        }
    }

    p[0] = 0; p[1] = nnz;
    i = (int*) realloc(i, nnz * sizeof(int));
    x = (T*)   realloc(x, nnz * sizeof(T));
    
    csc.m = m;
    csc.n = n;
    csc.nnz = nnz;
    csc.p = p;
    csc.i = i;
    csc.x = x;
}


template <typename T>
void coo_to_csc(
    const COO<T>& coo,
    CSC<T>& csc)
{
    int m = coo.m;
    int n = coo.n;
    int nnz = coo.nnz;

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



#endif