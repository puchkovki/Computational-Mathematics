// Implementation of Jacobi method for SLAE
#include "MatrixCSR.hpp"        // MatrixCSR class - compressed sparce row matrices
#include "Vector.hpp"           // Vector class - standard multidimension linear algebra vector
#include "Norms.hpp"            // Norm type
#include "Jacobi.hpp"           // Own header

// Makes first k iterations of Jacobi method
Vector JacobiK(const MatrixCSR& A, const Vector& b, size_t k) {

    Vector x_0 = b;
    Vector diag = A.get_Diag();
    size_t dim = A.get_Dim();

    for (size_t i = 0; i < dim; ++i) {
        x_0[i] /= diag[i];
    }

    Vector r = A * x_0 - b;
    MatrixCSR B = A.B_Jac();
    Vector d = x_0;

    for (size_t m = 0; m < k; ++m) {
        x_0 = B * x_0 + d;
        r = A * x_0 - b;
    }

    return x_0;
}

// Solves SLAE with accuracy eps by norm Norm, using Jacobi method
Vector Jacobi(const MatrixCSR& A, const Vector& b, norm Norm, double eps) {

    Vector x_0 = b;
    Vector diag = A.get_Diag();
    size_t dim = A.get_Dim();

    for (size_t i = 0; i < dim; ++i) {
        x_0[i] /= diag[i];
    }

    Vector r = A * x_0 - b;
    MatrixCSR B = A.B_Jac();
    Vector d = x_0;

    while (Norm(r) >= eps) {
        x_0 = B * x_0 + d;
        r = A * x_0 - b;
    }

    return x_0;
}

// Calculates the dependence of accuracy on iteration number for Jacobi method
Vector JacobiGraph(const MatrixCSR& A, const Vector& b, norm Norm, double eps) {

    Vector x_0 = b;
    Vector diag = A.get_Diag();
    size_t dim = A.get_Dim();

    for (size_t i = 0; i < dim; ++i) {
        x_0[i] /= diag[i];
    }

    Vector r = A * x_0 - b;
    MatrixCSR B = A.B_Jac();
    Vector d = x_0;
    Vector graph;

    while (Norm(r) >= eps) {
        graph.push_back(Norm(r));
        x_0 = B * x_0 + d;
        r = A * x_0 - b;
    }

    graph.push_back(Norm(r));
    
    return graph;
}