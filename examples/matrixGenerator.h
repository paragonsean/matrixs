#ifndef MATRIXGENERATOR_H
#define MATRIXGENERATOR_H

#include "../include/matrix.h"
#include "../include/matrix_storage_cep.h"
#include <string>
#include <iostream>
#include <fstream>

using namespace pnmatrix;

// Define a unique alias for the matrix type
typedef matrix<matrix_storage_cep<double>> MatrixType;

class MatrixGenerator {
public:
    static void generateRandomMatrixAndWriteToFile(const std::string& filename);
    static void readMatrixFromFile(MatrixType& A, MatrixType& b, const std::string& filename);
    static void writeOutSolution(const MatrixType& x, const std::string& filename);
    static void writeMatrixToFile(const MatrixType& A, const std::string& filename);
    static void printMatrix(const MatrixType& M, const std::string& name = "Matrix");
};

#endif // MATRIXGENERATOR_H
