#include "matrixGenerator.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <random>
#include <ctime>
#include "../include/matrix.h"
#include "../include/matrix_storage_cep.h"

typedef pnmatrix::matrix<pnmatrix::matrix_storage_cep<double>> matrix; // Replace with your actual container type
void MatrixGenerator::generateRandomMatrixAndWriteToFile(const std::string& filename) {
    int n;
    std::cout << "Enter the size of the matrix (n): ";
    while (!(std::cin >> n) || n <= 0) {
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout << "Invalid input, please enter a positive integer: ";
    }

    // Create matrices A (n x n), x (n x 1), and b (n x 1)
    pnmatrix::matrix<pnmatrix::matrix_storage_cep<double>> A(n, n);
    pnmatrix::matrix<pnmatrix::matrix_storage_cep<double>> x(n, 1);
    pnmatrix::matrix<pnmatrix::matrix_storage_cep<double>> b(n, 1);

    // Initialize random number generator
    std::default_random_engine gen(static_cast<unsigned>(time(nullptr)));
    std::normal_distribution<double> ndist(0.0, 50.0);

    // Initialize x with values from 1 to n (0-based indexing)
    for (int i = 0; i < n; ++i) {
        x.set_value(i, 0, static_cast<double>(i + 1));
    }

    // Generate random matrix A and ensure diagonal dominance
    for (int i = 0; i < n; ++i) {
        double s = 0.0;
        for (int j = 0; j < n; ++j) {
            double r = ndist(gen);
            s += std::abs(r);
            A.set_value(i, j, r);
        }
        // Ensure diagonal dominance by setting the diagonal element to the sum of abs values
        A.set_value(i, i, s);
    }

    // Compute b = A * x
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
            sum += A.get_value(i, j) * x.get_value(j, 0);
        }
        b.set_value(i, 0, sum);
    }

    // Write matrix A and vector b to the file as augmented matrix form
    std::ofstream fout(filename);
    if (!fout) {
        std::cerr << "Error opening file " << filename << " for writing.\n";
        return;
    }

    // Write dimensions of A only
    fout << n << " " << n << "\n";
    // Write A and b together on the same lines
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            fout << A.get_value(i, j) << " ";
        }
        // Append the corresponding element of b at the end of the row
        fout << b.get_value(i, 0) << "\n";
    }

    fout.close();

    std::cout << "\nMatrix A and vector b have been written to " << filename << " as an augmented matrix.\n";
}


void MatrixGenerator::readMatrixFromFile(MatrixType& A, MatrixType& b, const std::string& filename) {
    std::ifstream input(filename);
    if (!input) {
        std::cerr << "Failed to open the file: " << filename << "\n";
        throw std::runtime_error("File not found");
    }

    try {
        size_t rows, cols;
        input >> rows >> cols; // Read dimensions of A
        A = pnmatrix::matrix<pnmatrix::matrix_storage_cep<double>>(rows, cols);

        // Read A
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                double value;
                input >> value;
                A.set_value(i, j, value);
            }
        }

        // Read dimensions of b
        input >> rows >> cols;
        if (cols != 1) {
            throw std::invalid_argument("Vector b must be a column vector.");
        }

        b = pnmatrix::matrix<pnmatrix::matrix_storage_cep<double>>(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            double value;
            input >> value;
            b.set_value(i, 0, value);
        }

        std::cout << "Matrix A and vector b successfully loaded from " << filename << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Error reading matrix from file: " << e.what() << "\n";
        throw;
    }

    input.close();
}

void MatrixGenerator::writeOutSolution(const MatrixType& x, const std::string& filename) {
    std::ofstream output(filename);
    if (!output) {
        std::cerr << "Error opening file " << filename << " for writing.\n";
        return;
    }

    output << x; // Write the solution vector
    output.close();
    std::cout << "Solution written to " << filename << "\n";
}

void MatrixGenerator::writeMatrixToFile(const MatrixType& A, const std::string& filename) {
    std::ofstream output(filename);
    if (!output) {
        std::cerr << "Error opening file " << filename << " for writing.\n";
        return;
    }

    output << A; // Write the matrix
    output.close();
    std::cout << "Matrix written to " << filename << "\n";
}

void MatrixGenerator::printMatrix(const MatrixType& M, const std::string& name) {
    std::cout << name << " (" << M.get_row() << "x" << M.get_column() << "):\n";
    for (size_t i = 0; i < M.get_row(); ++i) {
        for (size_t j = 0; j < M.get_column(); ++j) {
            std::cout << std::setw(10) << std::setprecision(5) << M.get_value(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << std::endl;
}
