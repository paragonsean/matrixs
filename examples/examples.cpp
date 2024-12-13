#include "../include/matrix.h"
#include "../include/matrix_storage_cep.h"
#include "../include/gmres_solver.h"
#include "../include/gauss_seidel_solver.h"
#include "../include/jacobian_solver.h"
#include <iostream>

using namespace pnmatrix;

// Simple print function to display matrix values with zero-based indexing
template< typename MatrixType>
void print_matrix(const MatrixType& m) {
    for(auto row = m.begin(); row != m.end(); ++row) {
        for(auto col = row.begin(); col != row.end(); ++col) {
            std::cout << "(" << col.row_index() << ", " << col.column_index() << ") " << *col << " ";
        }
        std::cout << "\n";
    }
}

void gmres_example() {
    matrix<matrix_storage_cep<double>> m(3, 3);
    // Zero-based indexing
    m.set_value(0, 0, 1);
    m.set_value(0, 1, 1);
    m.set_value(0, 2, 1);
    m.set_value(1, 0, 0);
    m.set_value(1, 1, 4);
    m.set_value(1, 2, -1);
    m.set_value(2, 0, 2);
    m.set_value(2, 1, -2);
    m.set_value(2, 2, 1);
    
    std::cout << "Example GMRES.\n";
    std::cout << "Matrix A : \n";
    print_matrix(m);
    
    matrix<matrix_storage_cep<double>> b(3, 1);
    // Zero-based indexing
    b.set_value(0, 0, 6);
    b.set_value(1, 0, 5);
    b.set_value(2, 0, 1);

    std::cout << "Matrix b : \n";
    print_matrix(b);
    std::cout << "Use GMRES method to solve A * x = b : \n";
    
    gmres::option op;
    op.m = 3;      // Set m to 3 for a 3x3 matrix
    op.rm = 1e-5;  // Set rm to 1e-5
    std::cout << "Restart m : " << op.m << ", Error tolerance rm : " << op.rm << std::endl;
    
    gmres solver(op);
    auto result = solver.solve(m, b);
    
    std::cout << "Result : x \n";
    print_matrix(result);
    std::cout << "###\n";
}

void jacobian_example() {
    matrix<matrix_storage_cep<double>> m(3, 3);
    // Zero-based indexing
    m.set_value(0, 0, 8);
    m.set_value(0, 1, -3);
    m.set_value(0, 2, 2);
    m.set_value(1, 0, 4);
    m.set_value(1, 1, 11);
    m.set_value(1, 2, -1);
    m.set_value(2, 0, 6);
    m.set_value(2, 1, 3);
    m.set_value(2, 2, 12);
    
    std::cout << "Example Jacobian.\n";
    std::cout << "Matrix A : \n";
    print_matrix(m);
    
    matrix<matrix_storage_cep<double>> b(3, 1);
    // Zero-based indexing
    b.set_value(0, 0, 20);
    b.set_value(1, 0, 33);
    b.set_value(2, 0, 36);
    
    std::cout << "Matrix b : \n";
    print_matrix(b);
    std::cout << "Use Jacobian method to solve A * x = b : \n";
    
    jacobian::option op;
    op.rm = 1e-6;
    std::cout << "rm : " << op.rm << std::endl;
    
    jacobian solver(op);
    auto result = solver.solve(m, b);
    
    std::cout << "Result : x \n";
    print_matrix(result);
    std::cout << "###\n";
}

void gauss_seidel_example() {
    matrix<matrix_storage_cep<double>> m(3, 3);
    // Zero-based indexing
    m.set_value(0, 0, 8);
    m.set_value(0, 1, -3);
    m.set_value(0, 2, 2);
    m.set_value(1, 0, 4);
    m.set_value(1, 1, 11);
    m.set_value(1, 2, -1);
    m.set_value(2, 0, 6);
    m.set_value(2, 1, 3);
    m.set_value(2, 2, 12);
    
    std::cout << "Example Gauss-Seidel.\n";
    std::cout << "Matrix A : \n";
    print_matrix(m);
    
    matrix<matrix_storage_cep<double>> b(3, 1);
    // Zero-based indexing
    b.set_value(0, 0, 20);
    b.set_value(1, 0, 33);
    b.set_value(2, 0, 36);
    
    std::cout << "Matrix b : \n";
    print_matrix(b);
    std::cout << "Use Gauss-Seidel method to solve A * x = b : \n";
    
    gauss_seidel::option op;
    op.rm = 1e-6;
    std::cout << "rm : " << op.rm << std::endl;
    
    gauss_seidel solver(op);
    auto result = solver.solve(m, b);
    
    std::cout << "Result : x \n";
    print_matrix(result);
    std::cout << "###\n";
}

int main(int argc, char* argv[]) {
    gmres_example();
    jacobian_example();
    gauss_seidel_example();
    return 0;
}
