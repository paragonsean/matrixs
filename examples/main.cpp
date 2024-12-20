#include <iostream>
#include <string>
#include <chrono>
#include <random>
#include <thread>
#include <fstream>
#include <vector>
#include "matrix_generator.h"
#include "gauss_seidel_solver.h"
#include "../include/dense_matrix_storage.h"
#include "../include/sparse_matrix_storage.h"
#include "../include/matrix.h"
#include "../include/gaussian_elimination.h" // Include Gaussian Elimination header

using namespace pnmatrix;

// Alias types for dense and sparse matrices
using DenseMatrix = matrix<dense_matrix_storage<double>>;
using SparseMatrix = matrix<sparse_matrix_storage<double>>;

// Function to print a matrix
template <typename MatrixType>
void print_matrix(const MatrixType& m) {
    for (auto row = m.begin(); row != m.end(); ++row) {
        for (auto col = row.begin(); col != row.end(); ++col) {
            std::cout << "(" << col.row_index() << ", " << col.column_index() << ") " << *col << " ";
        }
        std::cout << "\n";
    }
}

// Multi-threaded row generator for DenseMatrix
void generateRows(DenseMatrix& A, const DenseMatrix& x, DenseMatrix& b,
                  int startRow, int endRow, unsigned threadSeed) {
    std::default_random_engine gen(threadSeed);
    std::normal_distribution<double> ndist(0.0, 50.0);

    int n = static_cast<int>(A.get_row());

    for (int i = startRow; i < endRow; ++i) {
        double s = 0.0;
        for (int j = 0; j < n; ++j) {
            double r = ndist(gen);
            s += std::abs(r);
            A.set_value(i, j, r);
        }
        A.set_value(i, i, s); // Ensure diagonal dominance

        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
            sum += A.get_value(i, j) * x.get_value(j, 0);
        }
        b.set_value(i, 0, sum);
    }
}

// Random matrix generator for DenseMatrix
void generateRandomMatrixAndWriteToFile(const std::string& filename) {
    int n;
    std::cout << "Enter the size of the matrix (n): ";
    while (!(std::cin >> n) || n <= 0) {
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout << "Invalid input, please enter a positive integer: ";
    }

    auto start = std::chrono::high_resolution_clock::now();

    DenseMatrix A(n, n);
    DenseMatrix x(n, 1);
    DenseMatrix b(n, 1);

    // Initialize vector x with 1 to n
    for (int i = 0; i < n; ++i) {
        x.set_value(i, 0, static_cast<double>(i + 1));
    }

    // Multi-threaded row generation
    unsigned numThreads = std::max(1u, std::thread::hardware_concurrency());
    std::vector<std::thread> threads;
    int rowsPerThread = n / numThreads;
    int remainder = n % numThreads;
    int startRow = 0;

    std::random_device rd;
    unsigned baseSeed = rd();

    for (unsigned t = 0; t < numThreads; ++t) {
        int endRow = startRow + rowsPerThread + (t < remainder ? 1 : 0);
        threads.emplace_back(generateRows, std::ref(A), std::ref(x), std::ref(b), startRow, endRow, baseSeed + t);
        startRow = endRow;
    }

    for (auto& th : threads) {
        th.join();
    }

    // Write matrix A and vector b to file in augmented format
    std::ofstream fout(filename);
    if (!fout) {
        std::cerr << "Error opening file " << filename << " for writing.\n";
        return;
    }

    fout << n << " " << n << "\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            fout << A.get_value(i, j) << " ";
        }
        fout << b.get_value(i, 0) << "\n";
    }

    fout.close();
    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Matrix and vector b written to " << filename << " in "
              << std::chrono::duration<double>(end - start).count() << " seconds.\n";
}

// Test matrix generation and solving system Ax = b
template <typename MatrixType>
void testMatrixGenerationAndSolve(const std::string& filename, const std::string& matrixType, bool printSolution) {
    try {
        MatrixType A, b;
        MatrixGenerator generator;

        auto totalStart = std::chrono::high_resolution_clock::now();

        std::cout << "\n--- " << matrixType << " Matrix Test ---\n";

        auto genStart = std::chrono::high_resolution_clock::now();
        generateRandomMatrixAndWriteToFile(filename);
        auto genEnd = std::chrono::high_resolution_clock::now();
        std::cout << matrixType << " Matrix Generation Time: "
                  << std::chrono::duration<double>(genEnd - genStart).count() << " seconds.\n";

        auto loadStart = std::chrono::high_resolution_clock::now();
        generator.loadDefaultAB(A, b, filename);
        auto loadEnd = std::chrono::high_resolution_clock::now();
        std::cout << "Matrix Load Time: "
                  << std::chrono::duration<double>(loadEnd - loadStart).count() << " seconds.\n";

        gauss_seidel::option op;
        op.rm = 1e-6;
        gauss_seidel solver(op);

        auto solveStart = std::chrono::high_resolution_clock::now();
        auto result = solver.solve(A, b);
        auto solveEnd = std::chrono::high_resolution_clock::now();
        std::cout << "Solver Execution Time: "
                  << std::chrono::duration<double>(solveEnd - solveStart).count() << " seconds.\n";

        auto totalEnd = std::chrono::high_resolution_clock::now();
        std::cout << "Total Time for " << matrixType << " Matrix: "
                  << std::chrono::duration<double>(totalEnd - totalStart).count() << " seconds.\n";

        if (printSolution) {
            std::cout << "\nSolution x:\n";
            print_matrix(result);
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }
}

// Function to test Gaussian Elimination
template <typename MatrixType>
void testGaussianElimination(const std::string& filename, const std::string& matrixType, bool printSolution) {
    try {
        MatrixType A, b;
        MatrixGenerator generator;

        auto totalStart = std::chrono::high_resolution_clock::now();

        std::cout << "\n--- " << matrixType << " Matrix Gaussian Elimination Test ---\n";

        auto genStart = std::chrono::high_resolution_clock::now();
        generateRandomMatrixAndWriteToFile(filename);
        auto genEnd = std::chrono::high_resolution_clock::now();
        std::cout << matrixType << " Matrix Generation Time: "
                  << std::chrono::duration<double>(genEnd - genStart).count() << " seconds.\n";

        auto loadStart = std::chrono::high_resolution_clock::now();
        generator.loadDefaultAB(A, b, filename);
        auto loadEnd = std::chrono::high_resolution_clock::now();
        std::cout << "Matrix Load Time: "
                  << std::chrono::duration<double>(loadEnd - loadStart).count() << " seconds.\n";

        gaussian_elimination::option op;
        op.rm = 1e-6;
        gaussian_elimination solver(op);

        auto solveStart = std::chrono::high_resolution_clock::now();
        auto result = solver.solve(A, b);
        auto solveEnd = std::chrono::high_resolution_clock::now();
        std::cout << "Solver Execution Time: "
                  << std::chrono::duration<double>(solveEnd - solveStart).count() << " seconds.\n";

        auto totalEnd = std::chrono::high_resolution_clock::now();
        std::cout << "Total Time for " << matrixType << " Matrix (Gaussian Elimination): "
                  << std::chrono::duration<double>(totalEnd - totalStart).count() << " seconds.\n";

        if (printSolution) {
            std::cout << "\nSolution x:\n";
            print_matrix(result);
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }
}

int main() {
    const std::string filename = "matdata.txt";
    bool printSolution = false;
    int solverChoice = 0;

    std::cout << "Select the solver to use:\n";
    std::cout << "1. Gauss-Seidel\n";
    std::cout << "2. Gaussian Elimination\n";
    std::cout << "Enter choice (1/2): ";
    std::cin >> solverChoice;

    std::cout << "Print solution? (1 = yes, 0 = no): ";
    std::cin >> printSolution;

    // Store the total time for all 30 runs
    double total_time = 0.0;

    if (solverChoice == 1) {
        std::cout << "\nStarting Dense and Sparse Matrix Tests with Gauss-Seidel...\n";

        for (int i = 0; i < 30; ++i) {
            auto start = std::chrono::high_resolution_clock::now();
            testMatrixGenerationAndSolve<DenseMatrix>(filename, "Dense", printSolution);
            auto end = std::chrono::high_resolution_clock::now();
            total_time += std::chrono::duration<double>(end - start).count();
            std::cout << "Run " << i + 1 << " time: "
                      << std::chrono::duration<double>(end - start).count() << " seconds\n";
        }
    } else if (solverChoice == 2) {
        std::cout << "\nStarting Dense and Sparse Matrix Tests with Gaussian Elimination...\n";

        for (int i = 0; i < 30; ++i) {
            auto start = std::chrono::high_resolution_clock::now();
            testGaussianElimination<DenseMatrix>(filename, "Dense", printSolution);
            auto end = std::chrono::high_resolution_clock::now();
            total_time += std::chrono::duration<double>(end - start).count();
            std::cout << "Run " << i + 1 << " time: "
                      << std::chrono::duration<double>(end - start).count() << " seconds\n";
        }
    } else {
        std::cerr << "Invalid choice. Exiting...\n";
    }

    std::cout << "\nAverage time for 30 runs: " << total_time / 30.0 << " seconds\n";
    std::cout << "\nAll tests completed.\n";
    return 0;
}
