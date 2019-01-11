#ifndef __EIGEN_H__
#define __EIGEN_H__


// Computes condition number ||A|| ||A^{-1}|| of a matrix
template <typename T>
int condition_number(const Eigen::SparseMatrix<T>& S) {
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(S);
    Eigen::SparseMatrix<double> I(S.rows(),S.cols());
    I.setIdentity();
    auto Sinv = solver.solve(I);
    return S.norm() * Sinv.norm();
}

// convert CSC format to Eigen Sparse
template <typename T>
void csc_to_eigensparse(const CSC<T>& csc, Eigen::SparseMatrix<T>& S) {
    vector<Triplet<double>> triplets;
    for (int j = 0; j < csc.m; ++j) {
        for (int p = csc.p[j]; p < csc.p[j+1]; ++p) {
            triplets.emplace_back(csc.i[p], j, csc.x[p]);
        }
    } 
    S.setFromTriplets(triplets.begin(), triplets.end());
}


#endif