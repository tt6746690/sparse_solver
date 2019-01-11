#ifndef __FORMATS_H__
#define __FORMATS_H__

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
    COO(const COO& coo);
    COO& operator=(const COO& coo);
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
private:
    void deallocate();
    void copy(const COO& from, COO& to);
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
    CSC(const CSC& csc);
    CSC& operator=(const CSC& csc);
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
private:
    void deallocate();
    void copy(const CSC& from, CSC& to);
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
#include <limits>

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
COO<T>::COO(const COO<T>& coo) {
    copy(coo, *this);
}

template <typename T>
COO<T>& COO<T>::operator=(const COO<T>& coo) {
    copy(coo, *this);
    return *this;
}

template <typename T>
void COO<T>::copy(const COO<T>& from, COO<T>& to) {
    auto nnz = from.nnz;
    to.m    =   from.m;
    to.n    =   from.n;
    to.nnz  =   from.nnz;
    deallocate();
    to.rowidx = (int*) malloc(nnz * sizeof(int));
    to.colidx = (int*) malloc(nnz * sizeof(int));
    to.values = (T*)   malloc(nnz * sizeof(T));
    std::copy(from.rowidx, from.rowidx + nnz, to.rowidx);
    std::copy(from.colidx, from.colidx + nnz, to.colidx);
    std::copy(from.values, from.values + nnz, to.values);
}

template <typename T>
void COO<T>::deallocate() {
    if (rowidx) free(rowidx);
    if (colidx) free(colidx);
    if (values) free(values);
}

template <typename T>
COO<T>::~COO() {
    deallocate();
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
CSC<T>::CSC(const CSC<T>& csc) {
    copy(csc, *this);
}

template <typename T>
CSC<T>& CSC<T>::operator=(const CSC<T>& csc) {
    copy(csc, *this);
    return *this;
}

template <typename T>
CSC<T>::~CSC()
{
#ifdef MEMORY_PRINT
    std::cout<<"CSC::~CSC()"<<": "<<this<<" ";
    if (!p && !i && !x) std::cout<<"ptrs are nullptr\n";
    else std::cout<<"releasing ptrs!\n";
#endif
    deallocate();
}

template <typename T>
void CSC<T>::copy(const CSC<T>& from, CSC<T>& to) {
    int m = from.m, n = from.n, nnz = from.nnz;
    to.m = m; to.n = n; to.nnz = nnz;
    deallocate();
    to.p = new int[n+1];
    to.i = new int[nnz];
    to.x = new T[nnz];
    std::copy(from.p, from.p + (n+1), to.p);
    std::copy(from.i, from.i + nnz,   to.i);
    std::copy(from.x, from.x + nnz,   to.x);
}

template <typename T>
void CSC<T>::deallocate() {
    if (p) delete[] p;
    if (i) delete[] i;
    if (x) delete[] x;
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
    size_t nnz = 0;

    for (int j = 0; j < v.size(); ++j) {
        if (v[j] != 0) {
            nnz += 1;
        }
    }

    auto p = new int[n+1];
    auto i = new int[nnz];
    auto x = new T[nnz];

    p[0] = 0; p[1] = nnz;

    nnz = 0;
    for (int j = 0; j < v.size(); ++j) {
        if (v[j] != 0) {
            if (nnz > std::numeric_limits<int>::max()) {
                printf("%zu (nnz) > %d (numerical limit for int)", nnz, std::numeric_limits<int>::max());
                exit(-1);
            }
            i[nnz] = j;
            x[nnz] = v[j];
            nnz += 1;
        }
    }
    
    csc.m = m;
    csc.n = n;
    csc.nnz = static_cast<int>(nnz);
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

    csc.p = new int[n+1];
    csc.i = new int[nnz];
    csc.x = new T[nnz];

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