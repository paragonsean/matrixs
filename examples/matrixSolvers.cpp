// matrixSolvers.cpp

#include "matrixSolvers.h"
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>        // for std::abs, std::max
#include <algorithm>    // for std::max
#include <thread>       // for std::thread
#include <mutex>        // for std::mutex
#include <atomic>       // for std::atomic
#include <functional>   // for std::ref
#include <sstream>      // for std::stringstream
#include "../include/matrix.h"
#include "../include/matrix_storage_cep.h"
#include "../include/gmres_solver.h"
#include "../include/gauss_seidel_solver.h"
#include "../include/jacobian_solver.h"

namespace pnmatrix {

typedef matrix<matrix_storage_cep<double>> MatrixType;

// Global tolerance for numerical precision
const double tolerance = 1e-6;

// Mutex for synchronizing output
std::mutex coutMutex;

// Function to safely print to std::cout
void safePrint(const std::string& message) {
    std::lock_guard<std::mutex> guard(coutMutex);
    std::cout << message;
}

// Function to compute relative error: ||Ax - b|| / ||b||
double computeRelativeError(const MatrixType& A, const MatrixType& x, const MatrixType& b) {
    int n = static_cast<int>(A.get_row());
    double numerator = 0.0;
    double denominator = 0.0;

    for (int i = 0; i < n; ++i) {
        double Ax_i = 0.0;
        for (int j = 0; j < static_cast<int>(A.get_column()); ++j) {
            Ax_i += A.get_value(i, j) * x.get_value(j, 0);
        }
        double residual = Ax_i - b.get_value(i, 0);
        numerator += residual * residual;
        denominator += b.get_value(i, 0) * b.get_value(i, 0);
    }

    return std::sqrt(numerator) / (std::sqrt(denominator) + tolerance);
}

int matrixSolvers::findPivot(const MatrixType& A, int column) {
    int pivotRow = column;
    double largestRatio = 0.0;

    for (int i = column; i < static_cast<int>(A.get_row()); ++i) {
        double rowMax = 0.0;
        for (size_t j = 0; j < A.get_column(); ++j) {
            rowMax = std::max(rowMax, std::abs(A.get_value(i, j)));
        }
        if (rowMax < tolerance) {
            continue;
        }

        double scaledPivot = std::abs(A.get_value(i, column)) / rowMax;
        if (scaledPivot > largestRatio) {
            largestRatio = scaledPivot;
            pivotRow = i;
        }
    }

    if (largestRatio < tolerance) {
        throw std::logic_error("Matrix is singular, cannot proceed with Gaussian elimination.");
    }

    return pivotRow;
}

void matrixSolvers::backsolve(const MatrixType& U, const MatrixType& y, MatrixType& x) {
    int n = static_cast<int>(U.get_row());

    if (static_cast<int>(x.get_row()) != n || x.get_column() != 1) {
        x = MatrixType(n, 1);
    }

    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += U.get_value(i, j) * x.get_value(j, 0);
        }
        double diag = U.get_value(i, i);
        if (std::abs(diag) < tolerance) {
            throw std::logic_error("Division by zero detected in backsolve.");
        }
        x.set_value(i, 0, (y.get_value(i, 0) - sum) / diag);
    }
}

void matrixSolvers::forwardSolve(MatrixType& L, MatrixType& y, const MatrixType& b) {
    int n = static_cast<int>(L.get_row());

    if (static_cast<int>(y.get_row()) != n || y.get_column() != 1) {
        y = MatrixType(n, 1);
    }

    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < i; ++j) {
            sum += L.get_value(i, j) * y.get_value(j, 0);
        }
        double diag = L.get_value(i, i);
        if (std::abs(diag) < tolerance) {
            throw std::logic_error("Division by zero detected in forwardSolve.");
        }
        y.set_value(i, 0, (b.get_value(i, 0) - sum) / diag);
    }
}

// Helper function for Gaussian Elimination row operations
void eliminateRowGauss(MatrixType& augmented, int pivotRow, int i, int j, 
                      std::atomic<int>& totalSteps, int printInterval) {
    double multiplier = augmented.get_value(j, pivotRow);
    int numColumns = augmented.get_column(); // Cache the number of columns

    for (int k = pivotRow; k < numColumns; ++k) {
        double updatedValue = augmented.get_value(j, k) - augmented.get_value(i, k) * multiplier;
        augmented.set_value(j, k, updatedValue);
    }

    // Increment the step count atomically
    int currentStep = totalSteps.fetch_add(1) + 1;

    // Check if it's time to print progress
    if (printInterval > 0 && currentStep % printInterval == 0) {
        // Get the current thread ID
        std::stringstream ss;
        ss << "Thread " << std::this_thread::get_id() 
           << " | Gaussian Elimination progress: " 
           << currentStep << " steps completed.\n";
        safePrint(ss.str());

        // Note: Computing relative error here is avoided for thread safety and performance
    }
}

// Helper function for LU Decomposition L computation
void computeLU(MatrixType& L, MatrixType& U, const MatrixType& A, int i, int k, 
              std::atomic<int>& luSteps, int printInterval) {
    double sum = 0.0;
    for (int m = 0; m < i; ++m) {
        sum += L.get_value(k, m) * U.get_value(m, i);
    }
    double value = (A.get_value(k, i) - sum) / U.get_value(i, i);
    L.set_value(k, i, value);

    // Increment the LU step count atomically
    int currentStep = luSteps.fetch_add(1) + 1;

    // Check if it's time to print progress
    if (printInterval > 0 && currentStep % printInterval == 0) {
        // Get the current thread ID
        std::stringstream ss;
        ss << "Thread " << std::this_thread::get_id() 
           << " | LU Decomposition progress: " 
           << currentStep << " steps completed.\n";
        safePrint(ss.str());

        // Note: Computing relative error here is avoided for thread safety and performance
    }
}

void matrixSolvers::gaussianElimination(MatrixType& A, MatrixType& x, const MatrixType& b) {
    if (A.get_column() != A.get_row() || A.get_row() != b.get_row()) {
        throw std::logic_error("Matrix size mismatch");
    }

    int n = static_cast<int>(A.get_row());
    MatrixType augmented(n, n + 1);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            augmented.set_value(i, j, A.get_value(i, j));
        }
        augmented.set_value(i, n, b.get_value(i, 0));
    }

    std::atomic<int> totalEliminationSteps(0);
    int printInterval = 1000; // Print progress every 100,000 steps

    // Initialize the next print step
    std::atomic<int> nextPrintStep(printInterval);

    for (int i = 0; i < n; ++i) {
        // Pivot selection and row swapping
        int pivot = findPivot(augmented, i);
        if (pivot != i) {
            augmented.element_row_transform_swap(i, pivot);
        }

        // Normalize the pivot row
        double pivotValue = augmented.get_value(i, i);
        if (std::abs(pivotValue) > tolerance) {
            augmented.element_row_transform_multi(i, 1.0 / pivotValue);
        }

        // Parallelize the elimination of lower rows using std::thread
        int startRow = i + 1;
        if (startRow >= n) {
            continue;
        }

        int rowsToUpdate = n - startRow;
        int numThreads = std::thread::hardware_concurrency();
        if (numThreads == 0) numThreads = 2; // Fallback to 2 threads if hardware_concurrency is not well-defined

        int rowsPerThread = rowsToUpdate / numThreads;
        int remainingRows = rowsToUpdate % numThreads;

        std::vector<std::thread> threads;
        int currentStart = startRow;

        for (int t = 0; t < numThreads; ++t) {
            int currentEnd = currentStart + rowsPerThread + (t < remainingRows ? 1 : 0);
            if (currentStart >= n) break;

            threads.emplace_back([&, currentStart, currentEnd]() {
                for (int j = currentStart; j < currentEnd; ++j) {
                    eliminateRowGauss(augmented, i, i, j, totalEliminationSteps, printInterval);
                }
            });

            currentStart = currentEnd;
        }

        for (auto& th : threads) {
            if (th.joinable()) {
                th.join();
            }
        }

        // Check if it's time to compute and print relative error
        while (totalEliminationSteps.load() >= nextPrintStep.load()) {
            int expected = nextPrintStep.load();
            if (nextPrintStep.compare_exchange_strong(expected, expected + printInterval)) {
                // Successfully updated, proceed to compute and print relative error

                // Extract current U and y from the augmented matrix
                MatrixType U(n, n);
                MatrixType y(n, 1);
                for (int row = 0; row < n; ++row) {
                    for (int col = 0; col < n; ++col) {
                        U.set_value(row, col, augmented.get_value(row, col));
                    }
                    y.set_value(row, 0, augmented.get_value(row, n));
                }

                // Solve Ux = y to get current x
                MatrixType currentX(n, 1);
                backsolve(U, y, currentX);

                // Compute relative error
                double relError = computeRelativeError(A, currentX, b);

                // Print relative error
                std::stringstream ss;
                ss << "Relative Error after " << totalEliminationSteps.load() 
                   << " steps: " << relError << "\n";
                safePrint(ss.str());
            } else {
                // Another thread has already handled this print step
                break;
            }
        }
    }

    // Extract U and y from the augmented matrix
    MatrixType U(n, n);
    MatrixType y(n, 1);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            U.set_value(i, j, augmented.get_value(i, j));
        }
        y.set_value(i, 0, augmented.get_value(i, n));
    }

    // Call backsolve with U and y
    backsolve(U, y, x);

    // Compute and print final relative error
    double finalRelError = computeRelativeError(A, x, b);
    safePrint("Final Relative Error in Gaussian Elimination: " + std::to_string(finalRelError) + "\n");
    safePrint("Total elimination steps in Gaussian Elimination: " + std::to_string(totalEliminationSteps.load()) + "\n");
}

void matrixSolvers::LUdecomposition(MatrixType& A, MatrixType& x, const MatrixType& b) {
    if (A.get_column() != A.get_row() || A.get_row() != b.get_row()) {
        throw std::logic_error("Matrix size mismatch");
    }

    int n = static_cast<int>(A.get_row());
    MatrixType L(n, n);
    MatrixType U(n, n);

    for (int i = 0; i < n; ++i) {
        L.set_value(i, i, 1.0);
    }

    std::atomic<int> luSteps(0);
    int printInterval = 100; // Print progress every 100,000 steps

    // Initialize the next print step
    std::atomic<int> nextPrintStep(printInterval);

    for (int i = 0; i < n; ++i) {
        // Compute U[i, i..n-1]
        for (int k = i; k < n; ++k) {
            luSteps.fetch_add(1);
            if (luSteps.load() >= nextPrintStep.load() && luSteps.load() % printInterval == 0) {
                // Attempt to update the nextPrintStep
                int expected = nextPrintStep.load();
                if (nextPrintStep.compare_exchange_strong(expected, expected + printInterval)) {
                    // Successfully updated, proceed to compute and print relative error

                    std::stringstream ss;
                    ss << "Main Thread | LU Decomposition progress: " 
                       << luSteps.load() << " steps completed.\n";
                    safePrint(ss.str());

                    // Forward substitution to solve Ly = b
                    MatrixType y(n, 1);
                    forwardSolve(L, y, b);

                    // Back substitution to solve Ux = y
                    MatrixType currentX(n, 1);
                    backsolve(U, y, currentX);

                    // Compute relative error
                    double relError = computeRelativeError(A, currentX, b);

                    // Print relative error
                    ss.str("");
                    ss << "Relative Error after " << luSteps.load() 
                       << " steps: " << relError << "\n";
                    safePrint(ss.str());
                }
            }

            double sum = 0.0;
            for (int m = 0; m < i; ++m) {
                sum += L.get_value(i, m) * U.get_value(m, k);
            }
            U.set_value(i, k, A.get_value(i, k) - sum);
        }

        double diag = U.get_value(i, i);
        if (std::abs(diag) < tolerance) {
            throw std::logic_error("LU decomposition failed: zero pivot encountered. Consider pivoting.");
        }

        // Compute L[i+1..n-1, i] in parallel using std::thread
        int startRow = i + 1;
        if (startRow < n) {
            int rowsToUpdate = n - startRow;
            int numThreads = std::thread::hardware_concurrency();
            if (numThreads == 0) numThreads = 2; // Fallback to 2 threads if hardware_concurrency is not well-defined

            int rowsPerThread = rowsToUpdate / numThreads;
            int remainingRows = rowsToUpdate % numThreads;

            std::vector<std::thread> threads;
            int currentStart = startRow;

            for (int t = 0; t < numThreads; ++t) {
                int currentEnd = currentStart + rowsPerThread + (t < remainingRows ? 1 : 0);
                if (currentStart >= n) break;

                threads.emplace_back([&, currentStart, currentEnd]() {
                    for (int k = currentStart; k < currentEnd; ++k) {
                        computeLU(L, U, A, i, k, luSteps, printInterval);
                    }
                });

                currentStart = currentEnd;
            }

            for (auto& th : threads) {
                if (th.joinable()) {
                    th.join();
                }
            }

            // Check if it's time to compute and print relative error
            while (luSteps.load() >= nextPrintStep.load()) {
                int expected = nextPrintStep.load();
                if (nextPrintStep.compare_exchange_strong(expected, expected + printInterval)) {
                    // Successfully updated, proceed to compute and print relative error

                    // Forward substitution to solve Ly = b
                    MatrixType y(n, 1);
                    forwardSolve(L, y, b);

                    // Back substitution to solve Ux = y
                    MatrixType currentX(n, 1);
                    backsolve(U, y, currentX);

                    // Compute relative error
                    double relError = computeRelativeError(A, currentX, b);

                    // Print relative error
                    std::stringstream ss;
                    ss << "Relative Error after " << luSteps.load() 
                       << " steps: " << relError << "\n";
                    safePrint(ss.str());
                } else {
                    // Another thread has already handled this print step
                    break;
                }
            }
        }
    }

    // Forward substitution to solve Ly = b
    MatrixType y(n, 1);
    forwardSolve(L, y, b);

    // Back substitution to solve Ux = y
    backsolve(U, y, x);

    // Compute and print final relative error
    double finalRelError = computeRelativeError(A, x, b);
    safePrint("Final Relative Error in LU Decomposition: " + std::to_string(finalRelError) + "\n");
    safePrint("Total LU decomposition steps: " + std::to_string(luSteps.load()) + "\n");
}

} // namespace pnmatrix
