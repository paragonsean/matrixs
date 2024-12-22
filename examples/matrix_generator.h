#ifndef MATRIX_GENERATOR_H
#define MATRIX_GENERATOR_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <random>
#include <ctime>
#include <thread>
#include <vector>
#include <atomic>
#include <sstream>
#include "../include/matrix.h"
#include "../include/dense_matrix_storage.h"
#include "../include/sparse_matrix_storage.h"

namespace pnmatrix {
class MatrixGenerator {
public:
    template <typename MatrixType>
    void generateRandomMatrixAndWriteToFile(const std::string& filename);

    template <typename MatrixType>
    void loadDefaultABTwo(MatrixType &A, MatrixType &b, const std::string &filename);

    template <typename MatrixType>
    void parallelParseLines(std::vector<std::string> &lines, MatrixType &A, MatrixType &b);


    template <typename MatrixType>
    void readMatrixFromFile(MatrixType &A, MatrixType &b, const std::string &filename);

    template <typename MatrixType>
    void print_matrix(const MatrixType &m);

    template <typename MatrixType>
    void loadDefaultAB(MatrixType &A, MatrixType &b, const std::string &filename);


    template <typename MatrixType>
    void writeMatrixToFile(const MatrixType& A, const std::string& filename);

    template <typename MatrixType>
    double calculate2NormError(const MatrixType& A, const MatrixType& x, const MatrixType& b);  
};



// Modified function to handle matrix output (added the writeMatrixToFile function)
template <typename MatrixType>
void MatrixGenerator::writeMatrixToFile(const MatrixType& A, const std::string& filename) {
    std::ofstream fout(filename);
    if (!fout) {
        std::cerr << "Error opening file " << filename << " for writing.\n";
        return;
    }

    fout << A.get_row() << " " << A.get_column() << "\n";
    for (size_t i = 0; i < A.get_row(); ++i) {
        for (size_t j = 0; j < A.get_column(); ++j) {
            fout << A.get_value(i, j) << " ";
        }
        fout << "\n";
    }
    fout.close();
}

template <typename MatrixType>
void MatrixGenerator::generateRandomMatrixAndWriteToFile(const std::string& filename) {
    int n;
    std::cout << "Enter the size of the matrix (n): ";
    while (!(std::cin >> n) || n <= 0) {
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout << "Invalid input, please enter a positive integer: ";
    }

    auto start = std::chrono::high_resolution_clock::now();
    MatrixType A(n, n);
    MatrixType x(n, 1);
    MatrixType b(n, 1);

    // Initialize x with values from 1 to n
    for (int i = 0; i < n; ++i) {
        x.set_value(i, 0, static_cast<double>(i + 1));
    }

    unsigned numThreads = std::min(std::max(1u, std::thread::hardware_concurrency()), static_cast<unsigned>(n));
    std::vector<std::thread> threads;
    int rowsPerThread = n / static_cast<int>(numThreads);
    int remainder = n % static_cast<int>(numThreads);
    int startRow = 0;

    std::random_device rd;
    unsigned baseSeed = rd();

    auto generateRows = [](MatrixType& A, const MatrixType& x, MatrixType& b, int startRow, int endRow, unsigned seed) {
        std::default_random_engine gen(seed);
        std::normal_distribution<double> ndist(0.0, 50.0);
        int n = static_cast<int>(A.get_row());
        for (int i = startRow; i < endRow; ++i) {
            double s = 0.0;
            for (int j = 0; j < n; ++j) {
                double r = ndist(gen);
                s += std::abs(r);
                A.set_value(i, j, r);
            }
            A.set_value(i, i, s);
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                sum += A.get_value(i, j) * x.get_value(j, 0);
            }
            b.set_value(i, 0, sum);
        }
    };

    for (unsigned t = 0; t < numThreads; ++t) {
        int endRow = startRow + rowsPerThread + ((static_cast<int>(t) < remainder) ? 1 : 0);
        if (startRow >= n) break;
        unsigned threadSeed = baseSeed + t;
        threads.emplace_back(generateRows, std::ref(A), std::ref(x), std::ref(b), startRow, endRow, threadSeed);
        startRow = endRow;
    }

    for (auto& th : threads) {
        th.join();
    }

    // Write matrix to file
    writeMatrixToFile(A, filename);

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Matrix A and vector b have been written to " << filename << " in augmented format.\n";
    std::cout << "Generation and writing took: "
              << std::chrono::duration<double>(end - start).count() << " seconds.\n";
}

// Templated function to parse lines into Matrix A and vector b
template <typename MatrixType>
void parseLines(const std::vector<std::string>& lines, MatrixType& A, MatrixType& b, size_t startRow, size_t endRow) {
    for (size_t i = startRow; i < endRow; ++i) {
        std::istringstream iss(lines[i]);
        for (size_t j = 0; j < A.get_column(); ++j) {
            double value;
            if (!(iss >> value)) {
                std::cerr << "Error: Invalid matrix value at row " << i << ", column " << j << "\n";
                return;
            }
            A.set_value(i, j, value);
        }

        double b_value;
        if (!(iss >> b_value)) {
            std::cerr << "Error: Invalid vector b value at row " << i << "\n";
            return;
        }
        b.set_value(i, 0, b_value);
    }
}

// Templated function to load Matrix A and vector b from a file
template <typename MatrixType>
void MatrixGenerator::loadDefaultABTwo(MatrixType& A, MatrixType& b, const std::string& filename) {
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

    A = MatrixType(rows, cols);
    b = MatrixType(rows, 1);

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

    // Parallel parsing of lines
    unsigned numThreads = std::max(1u, std::thread::hardware_concurrency());
    std::vector<std::thread> threads;
    size_t rowsPerThread = rows / numThreads;
    size_t remainder = rows % numThreads;
    size_t startRow = 0;

    for (unsigned t = 0; t < numThreads; ++t) {
        size_t endRow = startRow + rowsPerThread + (t < remainder ? 1 : 0);
        if (startRow >= rows) break;
        threads.emplace_back(parseLines<MatrixType>, std::ref(lines), std::ref(A), std::ref(b), startRow, endRow);
        startRow = endRow;
    }

    for (auto& th : threads) {
        th.join();
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Matrix A and vector b successfully loaded from " << filename << "\n";
    std::cout << "Loading took: " << std::chrono::duration<double>(end - start).count() << " seconds.\n";
}
// Parallel parsing of lines
template <typename MatrixType>
void MatrixGenerator::parallelParseLines(std::vector<std::string>& lines, MatrixType& A, MatrixType& b) {
    auto parseChunk = [&](size_t start, size_t end) {
        for (size_t i = start; i < end; ++i) {
            std::istringstream iss(lines[i]);
            for (size_t j = 0; j < A.get_column(); ++j) {
                double val;
                iss >> val;
                A.set_value(i, j, val);
            }
            double b_val;
            iss >> b_val;
            b.set_value(i, 0, b_val);
        }
    };

    unsigned numThreads = std::min(std::thread::hardware_concurrency(), static_cast<unsigned>(lines.size()));
    std::vector<std::thread> threads;
    size_t chunkSize = (lines.size() + numThreads - 1) / numThreads;

    for (unsigned t = 0; t < numThreads; ++t) {
        size_t start = t * chunkSize;
        size_t end = std::min(start + chunkSize, lines.size());
        threads.emplace_back(parseChunk, start, end);
    }

    for (auto& th : threads) th.join();
}
template <typename MatrixType>
void MatrixGenerator::readMatrixFromFile(MatrixType& A, MatrixType& b, const std::string& filename) {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        throw std::runtime_error("Error: Could not open file " + filename);
    }

    size_t rows, cols;
    inputFile >> rows >> cols;

    A = MatrixType(rows, cols);
    b = MatrixType(rows, 1);

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            double val;
            inputFile >> val;
            A.set_value(i, j, val);
        }
        double b_val;
        inputFile >> b_val;
        b.set_value(i, 0, b_val);
    }
}


// Simple print function to display matrix values with zero-based indexing
template< typename MatrixType>
void MatrixGenerator::print_matrix(const MatrixType& m) {
    for(auto row = m.begin(); row != m.end(); ++row) {
        for(auto col = row.begin(); col != row.end(); ++col) {
            std::cout << "(" << col.row_index() << ", " << col.column_index() << ") " << *col << " ";
        }
        std::cout << "\n";
    }
}
template <typename MatrixType>
void MatrixGenerator::loadDefaultAB(MatrixType& A, MatrixType& b, const std::string& filename) {
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

    A = MatrixType(rows, cols);
    b = MatrixType(rows, 1);
    std::vector<std::string> lines(rows);
    inputFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    for (size_t i = 0; i < rows; ++i) {
        if (!std::getline(inputFile, lines[i])) {
            std::cerr << "Error reading line " << i << " from file.\n";
            return;
        }
    }
    inputFile.close();
    parallelParseLines(lines, A, b);

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Matrix A and vector b successfully loaded from " << filename << "\n";
    std::cout << "Loading took: " << std::chrono::duration<double>(end - start).count() << " seconds.\n";
}
}
#endif // MATRIX_GENERATOR_H
