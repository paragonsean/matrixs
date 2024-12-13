// #include "Matrix<double> <double> .h"
// #include <iomanip>  // For std::setprecision and std::setw
// #include <algorithm>  // For std::max_element
// #include <vector>

// Matrix<double> <double> ::Matrix<double> <double> () {}

// Matrix<double> <double> ::Matrix<double> <double> (int rows, int cols) : array(rows) {
//     for (auto& row : array) {
//         row.resize(cols);
//     }
// }

// Matrix<double> <double> ::Matrix<double> <double> (int rows, int cols, double defaultValue) : array(rows) {
//     for (auto& row : array) {
//         row.assign(cols, defaultValue);
//     }
// }

// // Non-const version
// std::vector<double>& Matrix<double> <double> ::operator[](int row) {
//     if (row < 0 || row >= array.size()) {
//         std::cerr << "Error: Row index " << row << " out of range (size: " << array.size() << ")" << std::endl;
//         throw std::out_of_range("Matrix<double> <double> <double>index out of range");
//     }
//     return array[row];
// }

// // Const version
// const std::vector<double>& Matrix<double> <double> ::operator[](int row) const {
//     if (row < 0 || row >= array.size()) {
//         std::cerr << "Error: Row index " << row << " out of range (size: " << array.size() << ")" << std::endl;
//         throw std::out_of_range("Matrix<double> <double> <double>index out of range");
//     }
//     return array[row];
// }

// Matrix<double> <double> ::Matrix<double> <double> (const std::vector<std::vector<double>>& v) : array(v) {}

// Matrix<double> <double> ::Matrix<double> <double> (const std::vector<std::vector<double>>& v, double defaultValue) : array(v) {}

// Matrix<double> <double> ::Matrix<double> <double> (const std::vector<double>& v, double defaultValue) {
//     array.push_back(v);
// }

// Matrix<double> <double> ::Matrix<double> <double> (const Matrix<double> <double> & m) : array(m.array) {}

// Matrix<double> <double> ::~Matrix<double> <double> () {
//     array.clear();
// }

// bool Matrix<double> <double> ::operator==(const Matrix<double> <double> & m) const {
//     return array == m.array;
// }
// double Matrix<double> <double> ::rowMax(int row) const {
//         if (row < 0 || row >= numrows()) {
//             throw std::out_of_range("Row index out of range");
//         }
//         return *std::max_element(array[row].begin(), array[row].end(), [](double a, double b) {
//             return std::abs(a) < std::abs(b);
//         });
//     }
// bool Matrix<double> <double> ::operator!=(const Matrix<double> <double> & m) const {
//     return !(*this == m);
// }

// // Matrix<double> <double> & Matrix<double> <double> ::operator=(const Matrix<double> <double> & right) {
// //     if (this != &right) { // Correct self-assignment check
// //         array = right.array;
// //     }
// //     return *this;
// // }
// // Assignment Operator
// Matrix<double> <double> & Matrix<double> <double> ::operator=(const Matrix<double> <double> & other) {
//     if (this != &other) { // Prevent self-assignment
//         this->array = other.array; // Deep copy of std::vector elements
//     }
//     return *this;
// }
// Matrix<double> <double> & Matrix<double> <double> ::multiplyRow(int rowNum, double scaler) {
//     if (rowNum < 0 || rowNum >= array.size()) {
//         throw std::out_of_range("parameters given are outside of Matrix<double> <double> <double>boundaries");
//     }

//     std::vector<double>& row = array.at(rowNum);
//     for (int i = 0; i < row.size(); ++i) {
//         row[i] = scaler * row[i];
//     }
//     return *this;
// }

// Matrix<double> <double> & Matrix<double> <double> ::addRows(int rowOne, int rowTwo, int destinationRow) {
//     if (rowOne < 0 || rowOne >= array.size() || rowTwo < 0 || rowTwo >= array.size() || destinationRow < 0 || destinationRow >= array.size()) {
//         throw std::out_of_range("parameters given are outside of Matrix<double> <double> <double>boundaries");
//     }

//     std::vector<double>& row1 = array.at(rowOne);
//     std::vector<double>& row2 = array.at(rowTwo);
//     std::vector<double>& writeRow = array.at(destinationRow);

//     for (int i = 0; i < row1.size(); ++i) {
//         writeRow[i] = row1[i] + row2[i];
//     }
//     return *this;
// }

// Matrix<double> <double> & Matrix<double> <double> ::swapRows(std::vector<double>& swapRow, int rowOne) {
//     if (rowOne < 0 || rowOne >= array.size()) {
//         throw std::out_of_range("parameters given are outside of Matrix<double> <double> <double>boundaries");
//     }

//     array.at(rowOne).swap(swapRow);
//     return *this;
// }

// void Matrix<double> <double> ::eye() {
//     if (array.size() != array.at(0).size()) {
//         throw std::logic_error("Matrix<double> <double> <double>size mismatch");
//     }

//     int tempSize = array.size();
//     array.clear();
//     array.resize(tempSize);

//     for (int i = 0; i < array.size(); ++i) {
//         array.at(i).assign(tempSize, 0.0);
//         (*this)[i][i] = 1.0;
//     }
// }

// void Matrix<double> <double> ::eye(int rows, int cols) {
//     array.clear();
//     array.resize(rows);

//     for (int i = 0; i < array.size(); ++i) {
//         array.at(i).assign(cols, 0.0);
//         (*this)[i][i] = 1.0;
//     }
// }

// // Only keep one augment function to avoid duplicate errors
// Matrix<double> <double> & Matrix<double> <double> ::augment(const Matrix<double> <double> & right) {
//     if (array.size() != right.numrows() || right.numcols() != 1) {
//         throw std::logic_error("Matrix<double> <double> <double>size mismatch");
//     }

//     for (int i = 0; i < array.size(); ++i) {
//         array.at(i).push_back(right[i][0]);
//     }

//     return *this;
// }

// // Only keep one swapRows function to avoid duplicate errors
// Matrix<double> <double> & Matrix<double> <double> ::swapRows(int rowOne, int rowTwo) {
//     if (rowOne < 0 || rowOne >= array.size() || rowTwo < 0 || rowTwo >= array.size()) {
//         throw std::out_of_range("parameters given are outside of Matrix<double> <double> <double>boundaries");
//     }

//     array.at(rowOne).swap(array.at(rowTwo));
//     return *this;
// }

// Matrix<double> <double> <double>operator*(const Matrix<double> <double> & left, const Matrix<double> <double> & right) {
//     if (left.numcols() != right.numrows()) {
//         throw std::logic_error("Matrix<double> <double> <double>size mismatch");
//     }
//     Matrix<double> <double> <double>solution(left.numrows(), right.numcols());
//     for (int i = 0; i < left.numrows(); ++i) {
//         for (int j = 0; j < right.numcols(); ++j) {
//             double tempValue = 0.0;
//             for (int k = 0; k < right.numrows(); ++k) {
//                 tempValue += (left[i][k]) * (right[k][j]);
//             }
//             solution[i][j] = tempValue;
//         }
//     }
//     return solution;
// }

// std::ostream& operator<<(std::ostream& out, const Matrix<double> <double> & m) {
//     for (int j = 0; j < m.numrows(); ++j) {
//         for (int i = 0; i < m.numcols(); ++i) {
//             if (i > 0)
//                 out << ' ';
//             out << std::setprecision(18) << std::setw(21);
//             out << m[j][i];
//         }
//         out << "\n";
//     }
//     return out;
// }

// double relError(const Matrix<double> <double> & A, const Matrix<double> <double> & x, const Matrix<double> <double> & b) {
//     if (!(A.numrows() == b.numrows()) || !(x.numcols() == b.numcols())) {
//         throw std::logic_error("Matrix<double> <double> <double>size mismatch");
//     }
//     Matrix<double> <double> <double>est = A * x;
//     double diffSumOfSquares = 0.0;
//     double actualSumOfSquares = 0.0;

//     for (int i = 0; i < est.numrows(); ++i) {
//         for (int j = 0; j < est.numcols(); ++j) {
//             double diffSquared = pow((est[i][j] - b[i][j]), 2.0);
//             double actualSquared = pow(b[i][j], 2.0);

//             diffSumOfSquares += diffSquared;
//             actualSumOfSquares += actualSquared;
//         }
//     }

//     double err = (sqrt(diffSumOfSquares) / sqrt(actualSumOfSquares));
//     return err;
// }

// void compareMatrices(const Matrix<double> <double> & b, const Matrix<double> <double> & b_new) {
//     if (b.numrows() != b_new.numrows() || b.numcols() != b_new.numcols()) {
//         std::cerr << "Error: Matrix<double> <double> <double>size mismatch!" << std::endl;
//         return;
//     }
//     std::string me; 
//     double threshold = 1e-8;
//     std::cout << "Comparing matrices element-wise:\n";

//     for (int i = 0; i < b.numrows(); ++i) {
//         for (int j = 0; j < b.numcols(); ++j) {
//             double diff = std::abs(b[i][j] - b_new[i][j]);
//             double relError = std::abs(diff / (b[i][j] == 0.0 ? 1.0 : b[i][j]));  // To avoid division by zero

//             std::cout << "b[" << i << "][" << j << "] = " << b[i][j]
//                       << ", b_new[" << i << "][" << j << "] = " << b_new[i][j]
//                       << ", diff = " << diff
//                       << ", relative error = " << relError;

//             if (relError < threshold) {
//                 std::cout << " (Close enough)" << std::endl;
//             } else {
//                 std::cout << " (Not close enough)" << std::endl;
          
//             }
//         }
//     }
// }

// // Function to multiply Matrix<double> <double> <double>A by vector x and compute b_new
// void multiplyAx(const Matrix<double> <double> & A, const Matrix<double> <double> & x, const Matrix<double> <double> & b) {
//     if (A.numcols() != x.numrows()) {
//         std::cerr << "Error: The number of columns in A must match the number of rows in x!" << std::endl;
//         return;
//     }

//     // Compute b_new as A * x
//     Matrix<double> <double> <double>b_new(A.numrows(), 1);
//     for (int i = 0; i < A.numrows(); ++i) {
//         double sum = 0.0;
//         for (int j = 0; j < A.numcols(); ++j) {
//             sum += A[i][j] * x[j][0];  // Assuming x is a column vector
//         }
//         b_new[i][0] = sum;
//     }

//     std::cout << "Computed b_new (A * x):\n" << b_new << std::endl;

//     // Compare b and b_new
//     compareMatrices(b, b_new);
// }



// // Function to compute the absolute 2-norm error
// double calculate2NormError(const Matrix<double> <double> & A, const Matrix<double> <double> & x, const Matrix<double> <double> & b)
// {
//     if (b.numrows() != A.numrows() || b.numcols() != x.numcols()) {
//         throw std::logic_error("Matrix<double> <double> <double>size mismatch between A, x, and b.");
//     }

//     // Step 1: Calculate b2 = A * x
//     Matrix<double> <double> <double>b2 = A * x;

//     // Step 2: Calculate the error vector (b2 - b)
//     double sumSquares = 0.0;
//     for (int i = 0; i < b.numrows(); ++i)
//     {
//         for (int j = 0; j < b.numcols(); ++j)
//         {
//             double error = b2[i][j] - b[i][j];
//             sumSquares += std::pow(error, 2);
//         }
//     }

//     // Step 3: Calculate the square root of the sum (2-norm)
//     double norm2 = std::sqrt(sumSquares);

//     // Step 4: Divide by the number of elements (N)
//     int N = b.numrows() * b.numcols();
//     double errorNorm = norm2 / N;
    
//     return errorNorm;
// }
