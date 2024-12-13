#include <iostream>
#include <fstream>
#include <stdexcept>
#include <limits>
#include <iomanip>
#include <string>
#include <chrono>
#include <random>
#include <thread>
#include <vector>
#include <atomic>

#include "matrixGenerator.h"
#include "matrixSolvers.h"
#include "IterativeSolvers.h"
#include "../include/matrix.h"
#include "../include/matrix_storage_cep.h"
#include "../include/gmres_solver.h"       // Make sure these are included
#include "../include/jacobian_solver.h"    // 
#include "../include/gauss_seidel_solver.h"/
using namespace pnmatrix;

// Modified printMatrix to handle large matrices
template <typename MatrixType>
void printMatrix(const MatrixType& M, const std::string& name = "Matrix") {
    if (M.get_column() > 100) {
        std::cout << name << " is too large to print ("
                  << M.get_row() << "x" << M.get_column() << ").\n"
                  << "It will be written to 'large_matrix_output.txt' instead.\n";

        try {
            MatrixGenerator::writeMatrixToFile(M, "large_matrix_output.txt");
            std::cout << "Matrix written to large_matrix_output.txt\n";
        } catch (const std::exception &e) {
            std::cerr << "Error writing large matrix to file: " << e.what() << "\n";
        }
        return;
    }

    std::cout << name << " (" << M.get_row() << "x" << M.get_column() << "):\n";
    for (size_t i = 0; i < M.get_row(); ++i) {
        for (size_t j = 0; j < M.get_column(); ++j) {
            std::cout << std::setw(10) << std::setprecision(15) << M(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << std::endl;
}

void compareMatrices(const matrix<matrix_storage_cep<double>>& b, const matrix<matrix_storage_cep<double>>& b_new) {
    if (b.get_row() != b_new.get_row() || b.get_column() != b_new.get_column()) {
        std::cerr << "Error: Matrix size mismatch!" << std::endl;
        return;
    }

    double threshold = 1e-8;
    std::cout << "Comparing matrices element-wise:\n";

    for (size_t i = 0; i < b.get_row(); ++i) {
        for (size_t j = 0; j < b.get_column(); ++j) {
            double diff = std::abs(b(i, j) - b_new(i, j));
            double denom = (b(i, j) == 0.0 ? 1.0 : b(i, j));
            double relError = std::abs(diff / denom);

            std::cout << "b[" << i << "][" << j << "] = " << b(i, j)
                      << ", b_new[" << i << "][" << j << "] = " << b_new(i, j)
                      << ", diff = " << diff
                      << ", relative error = " << relError;

            if (relError < threshold) {
                std::cout << " (Close enough)" << std::endl;
            } else {
                std::cout << " (Not close enough)" << std::endl;
            }
        }
    }
}

void multiplyAx(const matrix<matrix_storage_cep<double>>& A, const matrix<matrix_storage_cep<double>>& x, const matrix<matrix_storage_cep<double>>& b) {
    if (A.get_column() != x.get_row()) {
        std::cerr << "Error: The number of columns in A must match the number of rows in x!" << std::endl;
        return;
    }

    matrix<matrix_storage_cep<double>> b_new(A.get_row(), 1);
    for (size_t i = 0; i < A.get_row(); ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < A.get_column(); ++j) {
            sum += A(i, j) * x(j, 0);
        }
        b_new.set_value(i, 0, sum);
    }

    std::cout << "Computed b_new (A * x):\n";
    printMatrix(b_new, "b_new");
    compareMatrices(b, b_new);
}

double calculate2NormError(const matrix<matrix_storage_cep<double>>& A, const matrix<matrix_storage_cep<double>>& x, const matrix<matrix_storage_cep<double>>& b)
{
    if (b.get_row() != A.get_row() || x.get_column() != b.get_column()) {
        throw std::logic_error("Size mismatch between A, x, and b.");
    }

    matrix<matrix_storage_cep<double>> b2 = A * x;
    double sumSquares = 0.0;
    for (size_t i = 0; i < b.get_row(); ++i) {
        for (size_t j = 0; j < b.get_column(); ++j) {
            double error = b2(i, j) - b(i, j);
            sumSquares += error * error;
        }
    }

    double norm2 = std::sqrt(sumSquares);
    int N = static_cast<int>(b.get_row() * b.get_column());
    double errorNorm = norm2 / N;

    return errorNorm;
}

void multiplyAndPrint(const matrix<matrix_storage_cep<double>>& A, const matrix<matrix_storage_cep<double>>& x) {
    if (A.get_column() != x.get_row()) {
        std::cerr << "Error: Matrix A and vector x dimensions are incompatible for multiplication.\n";
        return;
    }

    matrix<matrix_storage_cep<double>> result(A.get_row(), 1);
    for (size_t i = 0; i < A.get_row(); ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < A.get_column(); ++j) {
            sum += A(i, j) * x(j, 0);
        }
        result.set_value(i, 0, sum);
    }

    printMatrix(result, "A * x");
}

void saveMatrixToFile(const matrix<matrix_storage_cep<double>>& A, const std::string& filename) {
    try {
        MatrixGenerator::writeMatrixToFile(A, filename);
        std::cout << "Matrix A successfully saved to " << filename << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Error saving matrix A: " << e.what() << "\n";
    }
}

void generateRows(pnmatrix::matrix<pnmatrix::matrix_storage_cep<double>> &A,
                  const pnmatrix::matrix<pnmatrix::matrix_storage_cep<double>> &x,
                  pnmatrix::matrix<pnmatrix::matrix_storage_cep<double>> &b,
                  int startRow, int endRow, unsigned threadSeed) 
{
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
        // Set diagonal for dominance
        A.set_value(i, i, s);

        // Compute b[i,0]
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
            sum += A.get_value(i, j) * x.get_value(j, 0);
        }
        b.set_value(i, 0, sum);
    }
}


void generateRandomMatrixAndWriteToFile(const std::string& filename)
{
    int n;
    std::cout << "Enter the size of the matrix (n): ";
    while (!(std::cin >> n) || n <= 0)
    {
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout << "Invalid input, please enter a positive integer: ";
    }

    auto start = std::chrono::high_resolution_clock::now();

    pnmatrix::matrix<pnmatrix::matrix_storage_cep<double>> A(n, n);
    pnmatrix::matrix<pnmatrix::matrix_storage_cep<double>> x(n, 1);
    pnmatrix::matrix<pnmatrix::matrix_storage_cep<double>> b(n, 1);

    // Initialize x with values from 1 to n
    for (int i = 0; i < n; ++i) {
        x.set_value(i, 0, static_cast<double>(i + 1));
    }

    // Determine number of threads
    unsigned numThreads = std::max(1u, std::thread::hardware_concurrency());
    std::vector<std::thread> threads;
    int rowsPerThread = n / static_cast<int>(numThreads);
    int remainder = n % static_cast<int>(numThreads);
    int startRow = 0;

    // Use a random device to seed each thread differently
    std::random_device rd;
    unsigned baseSeed = rd();

    for (unsigned t = 0; t < numThreads; ++t) {
        int endRow = startRow + rowsPerThread + ((static_cast<int>(t) < remainder) ? 1 : 0);
        if (startRow >= n) break;
        unsigned threadSeed = baseSeed + t;
        threads.emplace_back(generateRows, std::ref(A), std::ref(x), std::ref(b),
                             startRow, endRow, threadSeed);
        startRow = endRow;
    }

    for (auto &th : threads) {
        th.join();
    }

    // Write A and b to file in augmented format
    std::ofstream fout(filename);
    if (!fout) {
        std::cerr << "Error opening file " << filename << " for writing.\n";
        return;
    }

    fout << n << " " << n << "\n";
    for (int i = 0; i < n; ++i)
    {
        for (int col = 0; col < n; ++col)
        {
            fout << A.get_value(i, col) << " ";
        }
        fout << b.get_value(i, 0) << "\n";
    }
    fout.close();

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "\nMatrix A and vector b have been written to " << filename << " in augmented format.\n";
    std::cout << "Generation and writing took: " 
              << std::chrono::duration<double>(end - start).count() << " seconds.\n";
}

// Parallel parsing of lines
static void parseLines(std::vector<std::string> &lines, matrix<matrix_storage_cep<double>> &A, matrix<matrix_storage_cep<double>> &b, size_t start, size_t end) {
    for (size_t i = start; i < end; ++i) {
        std::istringstream iss(lines[i]);
        for (size_t j = 0; j < A.get_column(); ++j) {
            double val;
            if (!(iss >> val)) {
                throw std::runtime_error("Failed to parse matrix A data.");
            }
            A.set_value(i, j, val);
        }
        double b_val;
        if (!(iss >> b_val)) {
            throw std::runtime_error("Failed to parse vector b data.");
        }
        b.set_value(i, 0, b_val);
    }
}

void loadDefaultAB(matrix<matrix_storage_cep<double>>& A, matrix<matrix_storage_cep<double>>& b, const std::string& filename) {
    auto start = std::chrono::high_resolution_clock::now();

    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << "\n";
        return;
    }

    size_t rows, cols;
    if (!(inputFile >> rows >> cols)) {
        std::cerr << "Failed to read dimensions for matrix A from file.\n";
        return;
    }
    if (rows == 0 || cols == 0) {
        std::cerr << "Matrix A dimensions must be greater than zero.\n";
        return;
    }

    A = matrix<matrix_storage_cep<double>>(rows, cols);
    b = matrix<matrix_storage_cep<double>>(rows, 1);

    // Read all lines first
    std::vector<std::string> lines(rows);
    inputFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // move to next line
    for (size_t i = 0; i < rows; ++i) {
        if (!std::getline(inputFile, lines[i])) {
            std::cerr << "Error reading line " << i << " from file.\n";
            return;
        }
    }

    inputFile.close();

    // Parallel parse
    unsigned numThreads = std::max(1u, std::thread::hardware_concurrency());
    std::vector<std::thread> threads;
    size_t rowsPerThread = rows / numThreads;
    size_t remainder = rows % numThreads;
    size_t startRow = 0;
    for (unsigned t = 0; t < numThreads; ++t) {
        size_t endRow = startRow + rowsPerThread + (t < remainder ? 1 : 0);
        if (startRow >= rows) break;
        threads.emplace_back(parseLines, std::ref(lines), std::ref(A), std::ref(b), startRow, endRow);
        startRow = endRow;
    }

    for (auto &th : threads) {
        th.join();
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Matrix A and vector b successfully loaded from " << filename << "\n";
    std::cout << "Loading took: " << std::chrono::duration<double>(end - start).count() << " seconds.\n";
}

int printSolverOptions(std::istream& input) {
    int selection;
    std::cout << "1 - Gaussian Elimination\n";
    std::cout << "2 - LU Decomposition\n";
    std::cout << "3 - Jacobi Iteration (From IterativeSolvers)\n";
    std::cout << "4 - Gauss-Seidel Iteration (From IterativeSolvers)\n";
    std::cout << "5 - SOR (From IterativeSolvers)\n";
    std::cout << "6 - GMRES (Custom Implementation)\n";
    std::cout << "7 - Jacobi (Custom Implementation)\n";
    std::cout << "8 - Gauss-Seidel (Custom Implementation)\n";
    std::cout << "Enter selection: ";
    while (!(input >> selection) || selection < 1 || selection > 8) {
        input.clear();
        input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout << "Invalid input, please try again: ";
    }
    return selection;
}

int executeSolver(int selection, matrix<matrix_storage_cep<double>>& A_original,  matrix<matrix_storage_cep<double>>& b, matrix<matrix_storage_cep<double>>& x_computed) {
    if (A_original.get_row() != A_original.get_column()) {
        std::cerr << "Error: Matrix A must be square.\n";
        return -1;
    }

    int N = static_cast<int>(A_original.get_row());
    x_computed = matrix<matrix_storage_cep<double>>(N, 1);

    double errLimit = 0.000001;
    int maxIterations = 100000;
    matrix<matrix_storage_cep<double>> A = A_original;

    auto t_start = std::chrono::high_resolution_clock::now();
    std::string methodName;

    try {
        switch (selection) {
            case 1:
                // Gaussian Elimination
                matrixSolvers msolver;
                msolver.gaussianElimination(A, x_computed, b);
                methodName = "Gaussian Elimination";
                break;
            case 2:
                // LU Decomposition
                matrixSolvers msolverlu;
                msolverlu.LUdecomposition(A, x_computed, b);
                methodName = "LU Decomposition";
                break;
            case 3:
                // Jacobi Iteration from IterativeSolvers
                std::cout << "Enter error limit (default 0.000001): ";
                std::cin >> errLimit;
                IterativeSolvers::jacobi(A, x_computed, b, errLimit, maxIterations);
                methodName = "Jacobi Iteration (IterativeSolvers)";
                break;
            case 4:
                // Gauss-Seidel Iteration from IterativeSolvers
                std::cout << "Enter error limit (default 0.000001): ";
                std::cin >> errLimit;
                IterativeSolvers::gaussSeidel(A, x_computed, b, errLimit, maxIterations);
                methodName = "Gauss-Seidel Iteration (IterativeSolvers)";
                break;
            case 5: {
                // SOR from IterativeSolvers
                double omega;
                std::cout << "Enter relaxation factor (omega): ";
                std::cin >> omega;
                std::cout << "Enter error limit (default 0.000001): ";
                std::cin >> errLimit;
                IterativeSolvers::sor(omega, A, x_computed, b, errLimit, maxIterations);
                methodName = "SOR (IterativeSolvers)";
                break;
            }
            case 6: {
                // GMRES (Custom Implementation)
                gmres::option op;
                std::cout << "Enter GMRES restart parameter m (default 30): ";
                std::cin >> op.m;
                std::cout << "Enter error limit rm (default 1e-6): ";
                std::cin >> op.rm;
                gmres solver(op);
                x_computed = solver.solve(A, b);
                methodName = "GMRES (Custom)";
                break;
            }
            case 7: {
                // Jacobi (Custom Implementation)
                jacobian::option op;
                std::cout << "Enter error limit rm (default 1e-6): ";
                std::cin >> op.rm;
                jacobian solver(op);
                x_computed = solver.solve(A, b);
                methodName = "Jacobi (Custom)";
                break;
            }
            case 8: {
                // Gauss-Seidel (Custom Implementation)
                gauss_seidel::option op;
                std::cout << "Enter error limit rm (default 1e-6): ";
                std::cin >> op.rm;
                gauss_seidel solver(op);
                x_computed = solver.solve(A, b);
                methodName = "Gauss-Seidel (Custom)";
                break;
            }
            default:
                return -1;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error during solving: " << e.what() << "\n";
        return -1;
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    std::cout << "Method: " << methodName << "\n";
    std::cout << "Calculation time: " << std::chrono::duration<double>(t_end - t_start).count() << " seconds.\n";
    if (x_computed.get_row() < 1000) {
        printMatrix(x_computed, "Solution");
    }

    return 0;
}

int main() {
    matrix<matrix_storage_cep<double>> A(1, 1), b(1, 1), x_computed(1, 1);
    int menuSelection = -1;

    while (menuSelection != 0) {
        std::cout << "\nMenu Options:\n";
        std::cout << "1 - Generate a random matrix A and vector b, and write to file (matdata.txt).\n";
        std::cout << "2 - Load matrix A and vector b from file (matdata.txt).\n";
        std::cout << "3 - Use a solver to solve for x.\n";
        std::cout << "4 - Print most recent computed solution x.\n";
        std::cout << "5 - Print matrix A.\n";
        std::cout << "6 - Print vector b.\n";
        std::cout << "7 - Multiply A * x and print the result.\n";
        std::cout << "8 - Write matrix A to file (matrix_A.txt).\n";
        std::cout << "0 - Exit.\n";
        std::cout << "Enter selection: ";

        while (!(std::cin >> menuSelection) || menuSelection < 0 || menuSelection > 8) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "Invalid input, please try again: ";
        }

        switch (menuSelection) {
            case 1:
                generateRandomMatrixAndWriteToFile("matdata.txt");
                break;
            case 2:
                loadDefaultAB(A, b, "matdata.txt");
                break;
            case 3:
                if (A.get_row() == 0 || b.get_row() == 0) {
                    std::cerr << "Error: Load matrix A and vector b first.\n";
                    break;
                }
                {
                    int solverSelection;
                    std::cout << "Select a solver:\n";
                    solverSelection = printSolverOptions(std::cin);

                    int status = executeSolver(solverSelection, A, b, x_computed);
                    if (status == 0 && x_computed.get_row() > 0 && x_computed.get_row() < 100) {
                        multiplyAx(A, x_computed, b);
                        double errorNorm = calculate2NormError(A, x_computed, b);
                        std::cout << "2-norm error after computation: " << errorNorm << "\n";
                    }
                }
                break;
            case 4:
                printMatrix(x_computed, "Solution x");
                break;
            case 5:
                printMatrix(A, "Matrix A");
                break;
            case 6:
                printMatrix(b, "Vector b");
                break;
            case 7:
                multiplyAndPrint(A, x_computed);
                break;
            case 8:
                saveMatrixToFile(A, "matrix_A.txt");
                break;
            case 0:
                std::cout << "Exiting program.\n";
                break;
            default:
                std::cerr << "Invalid selection.\n";
        }
    }

    return 0;
}
