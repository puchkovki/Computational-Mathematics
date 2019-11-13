// Definition of Jacobi method for SLAE
#pragma once

#include "MatrixCSR.hpp"        // MatrixCSR class - compressed sparce row matrices
#include "Vector.hpp"           // Vector class - standard multidimension linear algebra vector
#include "Norms.hpp"            // Norm type

// Makes first k iterations of Jacobi method
Vector JacobiK(const MatrixCSR& A, const Vector& b, size_t k);

// Solves SLAE with accuracy eps by norm Norm, using Jacobi method
Vector Jacobi(const MatrixCSR& A, const Vector& b, norm Norm, double eps);

// Calculates the dependence of accuracy on iteration number for Jacobi method
Vector JacobiGraph(const MatrixCSR& A, const Vector& b, norm Norm, double eps);