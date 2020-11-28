// Implementation of conjugated gradient method method for SLAE
#include "MatrixCSR.hpp"        // MatrixCSR class - compressed sparce row matrices
#include "Vector.hpp"           // Vector class - standard multidimension linear algebra vector
#include "Norms.hpp"            // Norm type
#include "ConjGradient.hpp"     // Own header

// Makes first k iterations of conjugated gradient method
Vector ConjGradientK(const MatrixCSR& A, const Vector& b, size_t k) {
    Vector x_0 = b;
    Vector r_0 = b - A * x_0;
    Vector r_1 = r_0;
    Vector z_0 = r_0;
    for (size_t m = 0; m < k; ++m) {
        double sq_0 = r_0 * r_0;
        Vector az = A * z_0;
        double alpha = sq_0 / (z_0 * az);
        x_0 = x_0  + alpha * z_0;
        r_1 = r_0 - alpha * az;
        double beta = (r_1 * r_1) / sq_0;
        z_0 = r_1 + beta * z_0;
        r_0 = r_1;
    }

    return x_0;
}

// Solves SLAE with accuracy eps by norm Norm, using conjugated gradient method
Vector ConjGradient(const MatrixCSR& A, const Vector& b, norm Norm, double eps) {
    Vector x_0 = b;
    Vector r_0 = b - A * x_0;
    Vector r_1 = r_0;
    Vector z_0 = r_0;
    while (Norm(r_0) >= eps) {
        double sq_0 = r_0 * r_0;
        Vector az = A * z_0;
        double alpha = sq_0 / (z_0 * az);
        x_0 = x_0  + alpha * z_0;
        r_1 = r_0 - alpha * az;
        double beta = (r_1 * r_1) / sq_0;
        z_0 = r_1 + beta * z_0;
        r_0 = r_1;
    }
    return x_0;
}

// Calculates the dependence of accuracy on iteration number for conjugated gradient method
Vector ConjGradientGraph(const MatrixCSR& A, const Vector& b, norm Norm, double eps) {
    Vector x_0 = b;
    Vector r_0 = b - A * x_0;
    Vector r_1 = r_0;
    Vector z_0 = r_0;
    Vector graph;
    while (Norm(r_0) >= eps) {
        graph.push_back(Norm(r_0));
        double sq_0 = r_0 * r_0;
        Vector az = A * z_0;
        double alpha = sq_0 / (z_0 * az);
        x_0 = x_0  + alpha * z_0;
        r_1 = r_0 - alpha * az;
        double beta = (r_1 * r_1) / sq_0;
        z_0 = r_1 + beta * z_0;
        r_0 = r_1;
    }
    graph.push_back(Norm(r_0));
    return graph;
}