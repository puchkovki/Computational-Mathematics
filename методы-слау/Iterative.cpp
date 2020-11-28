// Implementation of simple single-parameter iterative method method for SLAE
#include "General.hpp"          // General macroses
#include "MatrixCSR.hpp"        // MatrixCSR class - compressed sparce row matrices
#include "Vector.hpp"           // Vector class - standard multidimension linear algebra vector
#include "Norms.hpp"            // Norm type
#include "Iterative.hpp"        // Own header

// Makes first k iterations of simple single-parameter iterative method
Vector IterativeK(const MatrixCSR& A, const Vector& b, size_t k) {
    double tau = 1.9 / A.Spectral_Radius(EPS * 1e8);
    Vector x_0 = b;
    Vector r = A * x_0 - b;

    for (size_t m = 0; m < k; ++m) {
        x_0 = x_0 - tau * r;
        r = A * x_0 - b;
    }
    return x_0;
}

// Solves SLAE with accuracy eps by norm Norm, using simple single-parameter iterative method
Vector Iterative(const MatrixCSR& A, const Vector& b, norm Norm, double eps) {
    double tau = 1.9 / A.Spectral_Radius(eps * 1e8);
    Vector x_0 = b;
    Vector r = A * x_0 - b;

    while (Norm(r) >= eps) {
        x_0 = x_0 - tau * r;
        r = A * x_0 - b;
    }
    return x_0;
}

// Calculates the dependence of accuracy on iteration number for simple single-parameter iterative method
Vector IterativeGraph(const MatrixCSR& A, const Vector& b, norm Norm, double eps) {
    double tau = 1.9 / A.Spectral_Radius(eps * 1e8);
    Vector x_0 = b;
    Vector r = A * x_0 - b;
    Vector graph;
    while (Norm(r) >= eps) {
        graph.push_back(Norm(r));
        x_0 = x_0 - tau * r;
        r = A * x_0 - b;
    }
    graph.push_back(Norm(r));
    return graph;
}