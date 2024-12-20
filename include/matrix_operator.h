// #pragma once

// #include "matrix.h"
// #include "operator_proxy.h" // Ensure this includes necessary proxy class definitions
// #include <type_traits>
// #include <cassert>
// #include <functional>
// #include <iomanip>
// #include <iostream>
// #include <fstream>

// namespace pnmatrix {

// // Operator+ Overloads
// template<typename Proxy1, typename Proxy2,
//     typename std::enable_if<
//         both_real<
//             std::disjunction<is_op_type<Proxy1>, is_dense_matrix<Proxy1>>::value,
//             std::disjunction<is_op_type<Proxy2>, is_dense_matrix<Proxy2>>::value
//         >::value, int>::type = 0>
// auto operator+(const Proxy1& m1, const Proxy2& m2) -> op_add<Proxy1, Proxy2> {
//     assert(m1.get_row() == m2.get_row());
//     assert(m1.get_column() == m2.get_column());
//     return op_add<Proxy1, Proxy2>(m1, m2, m1.get_row(), m1.get_column());
// }

// template<typename MatrixType, typename std::enable_if<is_sparse_matrix<MatrixType>::value, int>::type = 0>
// auto operator+(const MatrixType& m1, const MatrixType& m2) -> MatrixType {
//     assert(m1.get_row() == m2.get_row() && m1.get_column() == m2.get_column());
//     MatrixType result(m1.get_row(), m1.get_column());
//     for (auto row_iter = m1.begin(); row_iter != m1.end(); ++row_iter) {
//         for (auto colu_iter = row_iter.begin(); colu_iter != row_iter.end(); ++colu_iter) {
//             result.set_value(colu_iter.row_index(), colu_iter.column_index(), *colu_iter);
//         }
//     }
//     for (auto row_iter = m2.begin(); row_iter != m2.end(); ++row_iter) {
//         for (auto colu_iter = row_iter.begin(); colu_iter != row_iter.end(); ++colu_iter) {
//             result.add_value(colu_iter.row_index(), colu_iter.column_index(), *colu_iter);
//         }
//     }
//     return result;
// }

// // Operator- Overloads
// template<typename Proxy1, typename Proxy2,
//     typename std::enable_if<
//         both_real<
//             std::disjunction<is_op_type<Proxy1>, is_dense_matrix<Proxy1>>::value,
//             std::disjunction<is_op_type<Proxy2>, is_dense_matrix<Proxy2>>::value
//         >::value, int>::type = 0>
// auto operator-(const Proxy1& m1, const Proxy2& m2) -> op_sub<Proxy1, Proxy2> {
//     assert(m1.get_row() == m2.get_row());
//     assert(m1.get_column() == m2.get_column());
//     return op_sub<Proxy1, Proxy2>(m1, m2, m1.get_row(), m1.get_column());
// }

// template<typename MatrixType, typename std::enable_if<is_sparse_matrix<MatrixType>::value, int>::type = 0>
// auto operator-(const MatrixType& m1, const MatrixType& m2) -> MatrixType {
//     assert(m1.get_row() == m2.get_row() && m1.get_column() == m2.get_column());
//     MatrixType result(m1.get_row(), m1.get_column());
//     for (auto row_iter = m1.begin(); row_iter != m1.end(); ++row_iter) {
//         for (auto colu_iter = row_iter.begin(); colu_iter != row_iter.end(); ++colu_iter) {
//             result.set_value(colu_iter.row_index(), colu_iter.column_index(), *colu_iter);
//         }
//     }
//     for (auto row_iter = m2.begin(); row_iter != m2.end(); ++row_iter) {
//         for (auto colu_iter = row_iter.begin(); colu_iter != row_iter.end(); ++colu_iter) {
//             result.add_value(colu_iter.row_index(), colu_iter.column_index(), -(*colu_iter));
//         }
//     }
//     return result;
// }

// // Operator* Overloads
// template<typename Proxy1, typename Proxy2,
//     typename std::enable_if<
//         both_real<
//             std::disjunction<is_op_type<Proxy1>, is_dense_matrix<Proxy1>>::value,
//             std::disjunction<is_op_type<Proxy2>, is_dense_matrix<Proxy2>>::value
//         >::value, int>::type = 0>
// auto operator*(const Proxy1& m1, const Proxy2& m2) -> op_mul<Proxy1, Proxy2> {
//     assert(m1.get_column() == m2.get_row());
//     return op_mul<Proxy1, Proxy2>(m1, m2, m1.get_row(), m2.get_column());
// }

// template<typename MatrixType, typename MatrixType2,
//     typename std::enable_if<
//         std::conjunction<is_sparse_matrix<MatrixType>, is_matrix_type<MatrixType2>>::value,
//     int>::type = 0>
// auto operator*(const MatrixType& m1, const MatrixType2& m2) -> MatrixType2 {
//     assert(m1.get_column() == m2.get_row());
//     using value_type = typename MatrixType::value_type;
//     static_assert(std::is_same<value_type, typename MatrixType2::value_type>::value, "error.");

//     MatrixType2 result(m1.get_row(), m2.get_column());
//     for (auto row = m1.begin(); row != m1.end(); ++row) {
//         for (size_type i = 0; i < m2.get_column(); ++i) {
//             value_type sum = value_type(0);
//             size_type row_ = row.row_index();
//             size_type colu_ = i;
//             for (auto col = row.begin(); col != row.end(); ++col) {
//                 sum += *col * m2.get_value(col.column_index(), i);
//             }
//             result.set_value(row_, colu_, sum);
//         }
//     }
//     return result;
// }

// // Scalar Division Overloads
// template<typename MatrixType, typename std::enable_if<is_sparse_matrix<MatrixType>::value, int>::type = 0>
// MatrixType operator/(const MatrixType& m, typename MatrixType::value_type value) {
//     MatrixType result(m.get_row(), m.get_column());
//     m.every_nozero_element([&](typename MatrixType::const_column_iterator iter) -> void {
//         result.set_value(iter.row_index(), iter.column_index(), *iter / value);
//     });
//     return result;
// }

// template<typename Proxy, typename std::enable_if<
//     std::disjunction<is_op_type<Proxy>, is_dense_matrix<Proxy>>::value, int>::type = 0>
// auto operator/(const Proxy& m, typename Proxy::value_type value) -> op_div_value<Proxy> {
//     return op_div_value<Proxy>(m, value, m.get_row(), m.get_column());
// }

// // Scalar Multiplication Overloads
// template<typename MatrixType, typename std::enable_if<is_sparse_matrix<MatrixType>::value, int>::type = 0>
// MatrixType operator*(const MatrixType& m, typename MatrixType::value_type value) {
//     MatrixType result(m.get_row(), m.get_column());
//     m.every_nozero_element([&](typename MatrixType::const_column_iterator iter) -> void {
//         result.set_value(iter.row_index(), iter.column_index(), *iter * value);
//     });
//     return result;
// }

// template<typename Proxy, typename std::enable_if<
//     std::disjunction<is_op_type<Proxy>, is_dense_matrix<Proxy>>::value, int>::type = 0>
// auto operator*(const Proxy& m, typename Proxy::value_type value) -> op_mul_value<Proxy> {
//     return op_mul_value<Proxy>(m, value, m.get_row(), m.get_column());
// }

// // Transpose Overloads
// template<typename Proxy, typename std::enable_if<
//     std::disjunction<is_op_type<Proxy>, is_dense_matrix<Proxy>>::value, int>::type = 0>
// op_tr<Proxy> tr(const Proxy& m) {
//     return op_tr<Proxy>(m, m.get_column(), m.get_row());
// }

// template<typename MatrixType, typename std::enable_if<is_sparse_matrix<MatrixType>::value, int>::type = 0>
// MatrixType tr(const MatrixType& m) {
//     MatrixType result(m.get_column(), m.get_row());
//     for (auto row_iter = m.begin(); row_iter != m.end(); ++row_iter) {
//         for (auto colu_iter = row_iter.begin(); colu_iter != row_iter.end(); ++colu_iter) {
//             result.set_value(colu_iter.column_index(), colu_iter.row_index(), *colu_iter);
//         }
//     }
//     return result;
// }

// // Operator<< Overload for Printing Matrices
// template <class Container>
// std::ostream& operator<<(std::ostream& out, const matrix<Container>& m) {
//     const size_t row_count = m.get_row();
//     const size_t column_count = m.get_column();

//     out << std::fixed << std::setprecision(4);

//     for (size_t i = 0; i < row_count; ++i) {
//         for (size_t j = 0; j < column_count; ++j) {
//             out << std::setw(10) << m.get_value(i, j);
//         }
//         out << '\n';
//     }
//     return out;
// }

// // Operator+ Overload Returning a Proxy (Avoid Duplication)
// template<typename Proxy1, typename Proxy2,
//     typename std::enable_if<
//         is_matrix_type<Proxy1>::value && is_matrix_type<Proxy2>::value &&
//         std::is_same<typename Proxy1::value_type, typename Proxy2::value_type>::value, int>::type = 0>
// auto operator+(const Proxy1& m1, const Proxy2& m2) -> op_add<Proxy1, Proxy2> {
//     if(m1.get_row() != m2.get_row() || m1.get_column() != m2.get_column()) {
//         throw std::invalid_argument("Matrix dimensions must agree for addition.");
//     }
//     return op_add<Proxy1, Proxy2>(m1, m2);
// }

// // Operator- Overload Returning a Proxy (Avoid Duplication)
// template<typename Proxy1, typename Proxy2,
//     typename std::enable_if<
//         is_matrix_type<Proxy1>::value && is_matrix_type<Proxy2>::value &&
//         std::is_same<typename Proxy1::value_type, typename Proxy2::value_type>::value, int>::type = 0>
// auto operator-(const Proxy1& m1, const Proxy2& m2) -> op_sub<Proxy1, Proxy2> {
//     if(m1.get_row() != m2.get_row() || m1.get_column() != m2.get_column()) {
//         throw std::invalid_argument("Matrix dimensions must agree for subtraction.");
//     }
//     return op_sub<Proxy1, Proxy2>(m1, m2);
// }

// // Operator* Overload Returning a Proxy for Matrix Multiplication (Avoid Duplication)
// template<typename Proxy1, typename Proxy2,
//     typename std::enable_if<
//         is_matrix_type<Proxy1>::value && is_matrix_type<Proxy2>::value &&
//         std::is_same<typename Proxy1::value_type, typename Proxy2::value_type>::value, int>::type = 0>
// auto operator*(const Proxy1& m1, const Proxy2& m2) -> op_mul<Proxy1, Proxy2> {
//     if(m1.get_column() != m2.get_row()) {
//         throw std::invalid_argument("Matrix dimensions must align for multiplication.");
//     }
//     return op_mul<Proxy1, Proxy2>(m1, m2);
// }

// } // namespace pnmatrix
