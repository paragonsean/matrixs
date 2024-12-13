

#ifndef MATRIXSOLVERS_H_
#define MATRIXSOLVERS_H_

#include "../include/matrix.h"
#include "../include/matrix_storage_cep.h"
#include "../include/gmres_solver.h"
#include "../include/gauss_seidel_solver.h"
#include "../include/jacobian_solver.h"

namespace pnmatrix {

typedef matrix<matrix_storage_cep<double>> MatrixType; // Define MatrixType using pnmatrix framework

class matrixSolvers {
public:
    /**
     * Finds the pivot row number in a given MatrixType.
     *
     * @param A        const reference to a diagonally dominant MatrixType.
     * @param column   column in which we are going to find the pivot row.
     * @return int     the row with the largest value to be used as the pivot.
     */
    int findPivot(const MatrixType& A, int column);
    void backsolve(const MatrixType& U, const MatrixType& y, MatrixType& x); // New declaration
    /**
     * Back substitution for a MatrixType in Upper Echelon form.
     *
     * @param A   reference to *augmented* MatrixType A, converted to REF or RREF.
     * @param x   reference to vector where the solution variables will be stored.
     */
    void backsolve(MatrixType& A, MatrixType& x);
    void backsolveUpper(const MatrixType& U, const MatrixType& y, MatrixType& x); // New declaration
    /**
     * Forward substitution for LU decomposition.
     *
     * @param L   reference to Lower triangular MatrixType where A = LU.
     * @param y   reference to intermediate vector y.
     * @param b   const reference to the RHS vector.
     */
    void forwardSolve(MatrixType& L, MatrixType& y, const MatrixType& b);

    /**
     * Performs Gaussian Elimination to solve a linear system.
     *
     * @param A   reference to square diagonally dominant MatrixType.
     * @param x   reference parameter where we will store the solution.
     * @param b   const reference - RHS, size A.numrows() x 1.
     */
    void gaussianElimination(MatrixType& A, MatrixType& x, const MatrixType& b);

    /**
     * Performs LU decomposition to solve a linear system.
     *
     * @param A   reference to square diagonally dominant MatrixType.
     * @param x   reference parameter where we will store the solution.
     * @param b   const reference - RHS, size A.numrows() x 1.
     */
    void LUdecomposition(MatrixType& A, MatrixType& x, const MatrixType& b);

    /**
     * Decomposes a square MatrixType into L and U components (in-place).
     *
     * @param A   reference to square MatrixType A to be decomposed.
     */
    void decomp(MatrixType& A);
};

} // namespace pnmatrix

#endif /* MATRIXSOLVERS_H_ */
